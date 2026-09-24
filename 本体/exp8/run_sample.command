#!/bin/zsh
RTD_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MPLCONFIGDIR="$RTD_SCRIPT_DIR/.matplotlib-cache" python3 "$RTD_SCRIPT_DIR/rtd.py" --time-unit s -o "$RTD_SCRIPT_DIR/sample_results" "$RTD_SCRIPT_DIR/examples/demo_three_tanks.csv"
if [[ $? -eq 0 ]]; then
  echo "示例处理完成。结果在 sample_results；按回车关闭。"
else
  echo "处理失败。请查看上面的错误信息；按回车关闭。"
fi
read
