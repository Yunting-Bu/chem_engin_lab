"""处理流动反应釜脉冲示踪实验的停留时间分布数据。

运行方式：python3 rtd.py 数据文件1.xls [数据文件2.xls ...]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _columns_from_headers(headers: list[str]) -> tuple[str, list[str]] | None:
    """按常见的时间、电导列名识别表格结构，保留探头在原表中的顺序。"""
    normalized = {
        name: re.sub(r"[\s()（）\[\]{}:：/_-]", "", name).lower()
        for name in headers
    }
    time_columns = [
        name for name, key in normalized.items()
        if re.fullmatch(r"(?:时间|time|t)(?:s|sec|秒|min|分钟|ms)?", key)
    ]
    signal_columns = [
        name for name, key in normalized.items()
        if key.startswith(("电导", "conductivity", "conductance", "κ"))
    ]
    if len(time_columns) == 1 and signal_columns:
        return time_columns[0], signal_columns
    return None


def _excel_column(ref: str) -> str:
    match = re.match(r"[A-Z]+", ref)
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {ref}")
    return match.group()


def _load_xlsx(path: Path, header_row: int | None) -> tuple[list[str], list[dict[str, Any]], int]:
    """读取第一张工作表，在连续数据表结束后停止。

    本次仪器的 .xls 文件实际是 OOXML 压缩包，并将格式延伸到第 65536 行。
    逐行读取可避免展开数百兆字节的空白格式单元格。
    """
    with zipfile.ZipFile(path) as archive:
        strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(item.itertext()) for item in strings_root]
        sheets = sorted(
            (name for name in archive.namelist()
             if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name)),
            key=lambda name: int(re.search(r"sheet(\d+)\.xml", name).group(1)),
        )
        if not sheets:
            raise ValueError(f"No worksheet found in {path}")

        headers: dict[str, str] = {}
        data: list[dict[str, Any]] = []
        blank_after_data = 0
        with archive.open(sheets[0]) as source:
            context = ET.iterparse(source, events=("start", "end"))
            _, root = next(context)
            for event, row in context:
                if event != "end" or row.tag != NS + "row":
                    continue
                number = int(row.get("r", "0"))
                values: dict[str, Any] = {}
                for cell in row.findall(NS + "c"):
                    col = _excel_column(cell.get("r", ""))
                    value_node = cell.find(NS + "v")
                    if value_node is not None and value_node.text is not None:
                        value: Any = value_node.text
                        if cell.get("t") == "s":
                            value = strings[int(value)]
                    else:
                        inline = cell.find(NS + "is")
                        value = "".join(inline.itertext()) if inline is not None else None
                    if value is not None and value != "":
                        values[col] = value

                if header_row is None and number <= 30:
                    candidate = {col: str(value).strip() for col, value in values.items()}
                    if _columns_from_headers(list(candidate.values())):
                        header_row = number
                if number == header_row:
                    headers = {col: str(value).strip() for col, value in values.items()}
                elif header_row is not None and number > header_row and headers:
                    record = {name: values.get(col) for col, name in headers.items()}
                    if any(value is not None for value in record.values()):
                        data.append(record)
                        blank_after_data = 0
                    elif data:
                        blank_after_data += 1
                        if blank_after_data >= 20:
                            break
                elif header_row is None and number > 30:
                    break
                row.clear()
                root.clear()

    if not headers or not data:
        raise ValueError(f"No contiguous data rows found in {path}")
    return list(headers.values()), data, header_row


def _load_csv(path: Path, header_row: int | None) -> tuple[list[str], list[dict[str, Any]], int]:
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            with path.open("r", encoding=encoding, newline="") as stream:
                delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
                all_rows = list(csv.reader(stream, delimiter=delimiter))
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError(f"Cannot decode CSV: {path}")
    if header_row is None:
        header_row = next(
            (index for index, cells in enumerate(all_rows[:30], start=1)
             if _columns_from_headers([cell.strip() for cell in cells])),
            None,
        )
    if header_row is None:
        raise ValueError(f"未在前 30 行找到时间列和电导列：{path}")
    if len(all_rows) <= header_row:
        raise ValueError(f"CSV has no data after header row {header_row}: {path}")
    headers = [cell.strip() for cell in all_rows[header_row - 1]]
    rows = []
    for cells in all_rows[header_row:]:
        if not any(cells):
            continue
        rows.append(dict(zip(headers, cells)))
    return headers, rows, header_row


def load_table(path: Path, header_row: int | None = None) -> tuple[list[str], list[dict[str, Any]], int]:
    if zipfile.is_zipfile(path):
        return _load_xlsx(path, header_row)
    if path.suffix.lower() in (".csv", ".tsv"):
        return _load_csv(path, header_row)
    raise ValueError(
        f"Unsupported file: {path}. Export genuine binary .xls files as .xlsx or CSV first."
    )


def _number(value: Any, description: str) -> float:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        raise ValueError(f"{description} is not numeric: {value!r}") from None
    if not math.isfinite(number):
        raise ValueError(f"{description} is not finite: {value!r}")
    return number


def detect_response_start(times: list[float], readings: list[float]) -> tuple[int, int]:
    """寻找第一探头主峰前最近的稳定水电导区间。

    从最大读数向前搜索，以免把实验开始时的仪器启动变化误判为示踪剂响应。
    返回的是最后一个稳定采样点，而不是未记录的实际注射瞬间。
    """
    if len(times) != len(readings) or len(times) < 40:
        raise ValueError("Automatic response detection needs at least 40 samples")
    peak_index = max(range(len(readings)), key=readings.__getitem__)
    width = min(30, max(10, peak_index // 3))
    if peak_index < width + 3:
        raise ValueError("No sufficiently long water-only segment before the main peak")
    for next_index in range(peak_index - 2, width, -1):
        window = readings[next_index - width:next_index]
        baseline = statistics.median(window)
        amplitude = readings[peak_index] - baseline
        step_noise = statistics.median(
            abs(b - a) for a, b in zip(window, window[1:])
        )
        tolerance = max(0.015 * amplitude, 6 * step_noise, abs(readings[peak_index]) * 1e-9, 1e-12)
        if (
            amplitude > 5 * tolerance
            and max(window) - min(window) <= tolerance
            and readings[next_index] > baseline + tolerance
            and max(readings[next_index:next_index + 3]) > baseline + 3 * tolerance
        ):
            return next_index - 1, width
    raise ValueError("Could not identify a stable water segment before a sustained response")


@dataclass
class Processed:
    label: str
    source: str
    sensor: str
    response_onset_time: float
    samples: list[dict[str, float]]
    summary: dict[str, float | str | None]


def process_signal(
    *, label: str, source: str, sensor: str, times: list[float],
    readings: list[float], start_index: int, baseline_width: int,
) -> Processed:
    if len(times) != len(readings) or len(times) < 5:
        raise ValueError("Need at least five time/conductivity pairs")
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("Time must be strictly increasing with no duplicates")
    if start_index < baseline_width - 1 or baseline_width < 3 or start_index >= len(times) - 4:
        raise ValueError("Invalid automatically detected response or baseline window")
    baseline_first = start_index - baseline_width + 1
    baseline_readings = readings[baseline_first:start_index + 1]
    baseline = float(statistics.median(baseline_readings))
    response_onset_time = times[start_index]
    selected = [(t - response_onset_time, t, y) for t, y in zip(times[start_index:], readings[start_index:])]
    if len(selected) < 5:
        raise ValueError("Too few points at or after the detected response onset")
    relative = [item[0] for item in selected]
    spacings = [b - a for a, b in zip(relative, relative[1:])]
    dt = statistics.median(spacings)
    if any(abs(gap - dt) > max(1e-6, abs(dt) * 1e-4) for gap in spacings):
        raise ValueError("Lecture discrete formula requires equally spaced time points")

    corrected = [max(y - baseline, 0.0) for _, _, y in selected]
    if max(corrected) <= 0:
        raise ValueError("All baseline-corrected readings are zero")
    total = sum(corrected) * dt
    first_moment = sum(t * c for t, c in zip(relative, corrected)) * dt
    second_moment = sum(t * t * c for t, c in zip(relative, corrected)) * dt
    mean = first_moment / total
    variance = second_moment / total - mean * mean
    if variance <= 0:
        raise ValueError("Nonpositive residence-time variance")
    effective_tanks = mean * mean / variance
    density = [c / total for c in corrected]
    cumulative = [0.0]
    for i in range(1, len(relative)):
        cumulative.append(
            cumulative[-1]
            + (density[i - 1] + density[i]) * (relative[i] - relative[i - 1]) / 2
        )
    peak = max(corrected)
    samples = [
        {
            "recorded_time": absolute,
            "t": t,
            "raw_conductivity": reading,
            "baseline": baseline,
            "c_proxy": c,
            "c_over_cmax": c / peak,
            "E_per_time": e,
            "F": f,
        }
        for (t, absolute, reading), c, e, f in zip(
            selected, corrected, density, cumulative
        )
    ]
    summary: dict[str, float | str | None] = {
        "label": label,
        "source": source,
        "sensor": sensor,
        "response_onset_time": response_onset_time,
        "time_zero_method": "last stable sample before first-probe rise",
        "water_baseline": baseline,
        "baseline_window_start": times[baseline_first],
        "baseline_window_end": times[start_index],
        "baseline_window_range": max(baseline_readings) - min(baseline_readings),
        "sample_interval": dt,
        "n_points": len(samples),
        "sum_c_dt": total,
        "sum_t_c_dt": first_moment,
        "sum_t2_c_dt": second_moment,
        "mean_time": mean,
        "variance": variance,
        "N": effective_tanks,
        "peak_c_proxy": peak,
        "final_c_proxy": corrected[-1],
        "final_fraction_of_peak": corrected[-1] / peak,
        "final_F": cumulative[-1],
    }
    return Processed(label, source, sensor, response_onset_time, samples, summary)


def _safe_filename(value: str) -> str:
    return re.sub(r"[^\w.-]+", "_", value).strip("_") or "run"


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _make_plots(series: list[Processed], output_dir: Path, time_unit: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.font_manager import FontProperties, fontManager
    except ImportError as exc:
        raise RuntimeError("Plotting needs matplotlib") from exc
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    chinese_font = Path("/System/Library/Fonts/STHeiti Medium.ttc")
    if not chinese_font.exists():
        candidates = (
            "Noto Sans CJK SC", "Noto Sans CJK JP", "Microsoft YaHei",
            "SimHei", "PingFang SC", "WenQuanYi Zen Hei", "Arial Unicode MS",
        )
        available = {font.name: Path(font.fname) for font in fontManager.ttflist}
        chinese_font = next((available[name] for name in candidates if name in available), None)
    legend_font = FontProperties(fname=str(chinese_font), size=8) if chinese_font else None
    if legend_font is None and any(
        not f"{result.label}: {result.sensor}".isascii() for result in series
    ):
        print("Warning: no Chinese plot font found; legends use Series 1, Series 2, etc.")
    charts = [
        ("c_over_cmax", "c(t)/c_max (relative concentration)", "c(t)-t response curve", "c_t.png"),
        ("E_per_time", f"E(t) (1/{time_unit})", "E(t)-t residence-time density", "E_t.png"),
        ("F", "F(t)", "F(t)-t cumulative distribution", "F_t.png"),
    ]
    for field, ylabel, title, filename in charts:
        fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
        for index, result in enumerate(series, start=1):
            proposed_label = f"{result.label}: {result.sensor}"
            plot_label = proposed_label if legend_font or proposed_label.isascii() else f"Series {index}"
            ax.plot(
                [row["t"] for row in result.samples],
                [row[field] for row in result.samples],
                label=plot_label, linewidth=1.8,
            )
        ax.set(title=title, xlabel=f"t ({time_unit}; since first response)", ylabel=ylabel)
        if field in ("c_over_cmax", "F"):
            ax.set_ylim(0, 1.04)
        ax.grid(alpha=0.2)
        ax.legend(prop=legend_font, fontsize=8 if legend_font is None else None, ncol=2)
        fig.savefig(output_dir / filename, dpi=180)
        plt.close(fig)


def _run_config(config: dict[str, Any], base_dir: Path) -> Path:
    if not isinstance(config.get("runs"), list) or not config["runs"]:
        raise ValueError("Config must contain a nonempty 'runs' list")
    output_dir = Path(config.get("output_dir", "rtd_output")).expanduser()
    if not output_dir.is_absolute():
        output_dir = base_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    processed: list[Processed] = []

    for run_config in config["runs"]:
        source = Path(run_config["file"]).expanduser()
        if not source.is_absolute():
            source = base_dir / source
        if not source.exists():
            raise FileNotFoundError(source)
        if source.resolve().parent == output_dir.resolve():
            raise ValueError("Input file and output directory must be separate to protect the source data")
        configured_header = run_config.get("header_row")
        header_row = int(configured_header) if configured_header is not None else None
        if header_row is not None and header_row < 1:
            raise ValueError("header_row must be at least 1")
        headers, rows, header_row = load_table(source, header_row)
        inferred_columns = _columns_from_headers(headers)
        time_col = run_config.get("time_column") or (
            inferred_columns[0] if inferred_columns else "时间"
        )
        if time_col not in headers:
            raise ValueError(f"Time column {time_col!r} not found in {source}; found {headers}")
        signal_cols = run_config.get("signal_columns")
        if signal_cols is None:
            signal_cols = inferred_columns[1] if inferred_columns else []
        if (not isinstance(signal_cols, list) or not signal_cols
                or any(not isinstance(col, str) or col not in headers for col in signal_cols)):
            raise ValueError(f"Signal columns not found in {source}; found {headers}")
        if len(set(signal_cols)) != len(signal_cols):
            raise ValueError("signal_columns contains duplicate names")
        times: list[float] = []
        signals: dict[str, list[float]] = {col: [] for col in signal_cols}
        for row_number, row in enumerate(rows, start=header_row + 1):
            if row.get(time_col) in (None, ""):
                continue
            times.append(_number(row[time_col], f"{source.name} row {row_number} time"))
            for col in signal_cols:
                signals[col].append(_number(row.get(col), f"{source.name} row {row_number} {col}"))

        label = run_config.get("label", source.stem)
        start_index, baseline_width = detect_response_start(times, signals[signal_cols[0]])
        for col in signal_cols:
            result = process_signal(
                label=label, source=str(source), sensor=col,
                times=times, readings=signals[col],
                start_index=start_index, baseline_width=baseline_width,
            )
            result.summary["time_unit"] = str(config.get("time_unit", "recorded unit"))
            processed.append(result)

    filenames = [_safe_filename(f"{p.label}_{p.sensor}") + ".csv" for p in processed]
    if len(set(filenames)) != len(filenames):
        raise ValueError("Two signals would create the same output filename; give runs distinct labels")
    for result, filename in zip(processed, filenames):
        _write_csv(output_dir / filename, result.samples)
    _write_csv(output_dir / "results.csv", [p.summary for p in processed])
    (output_dir / "results.json").write_text(
        json.dumps([p.summary for p in processed], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    _make_plots(processed, output_dir, str(config.get("time_unit", "recorded unit")))
    return output_dir


def run(config_path: Path) -> Path:
    """兼容已有的 JSON 配置文件；普通用户可直接传数据文件。"""
    config = json.loads(config_path.read_text(encoding="utf-8"))
    return _run_config(config, config_path.parent)


def run_files(paths: list[Path], output_dir: Path | None = None, time_unit: str = "recorded unit") -> Path:
    """直接从一个或多个实验数据文件推断表头、列名和水电导。"""
    sources = [path.expanduser().resolve() for path in paths]
    if not sources:
        raise ValueError("请至少提供一个实验数据文件")
    destination = (output_dir.expanduser().resolve() if output_dir
                   else sources[0].parent / "rtd_results")
    counts: dict[str, int] = {}
    runs = []
    for source in sources:
        counts[source.stem] = counts.get(source.stem, 0) + 1
        suffix = f"_{counts[source.stem]}" if counts[source.stem] > 1 else ""
        runs.append({"file": str(source), "label": source.stem + suffix})
    return _run_config({
        "output_dir": str(destination), "time_unit": time_unit, "runs": runs,
    }, Path.cwd())


def main() -> int:
    parser = argparse.ArgumentParser(description="处理实验 8 的脉冲示踪电导数据")
    parser.add_argument("inputs", nargs="+", type=Path, help="一个或多个 CSV、TSV、XLSX 或 XLS 数据文件")
    parser.add_argument("-o", "--output", type=Path, help="结果目录；默认在第一份数据旁创建 rtd_results")
    parser.add_argument("--time-unit", default="recorded unit", help="图中时间单位，例如 s；默认使用记录单位")
    args = parser.parse_args()
    try:
        if len(args.inputs) == 1 and args.inputs[0].suffix.lower() == ".json":
            output_dir = run(args.inputs[0].expanduser().resolve())
        elif any(path.suffix.lower() == ".json" for path in args.inputs):
            raise ValueError("JSON 配置文件不能与数据文件同时传入")
        else:
            output_dir = run_files(args.inputs, args.output, args.time_unit)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"Results and plots saved to: {output_dir}")
    with (output_dir / "results.json").open(encoding="utf-8") as stream:
        summaries = json.load(stream)
    for item in summaries:
        print(
            f"{item['label']} {item['sensor']}: response onset={item['response_onset_time']:.3f}, "
            f"water baseline={item['water_baseline']:.3f}, "
            f"mean={item['mean_time']:.3f}, variance={item['variance']:.3f}, N={item['N']:.3f}"
        )
        if item["final_fraction_of_peak"] > 0.05:
            print("  Warning: the measured curve still has a substantial tail at its final point.")
        if item["N"] < 1:
            print("  Warning: N < 1 is outside the ideal tanks-in-series interpretation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
