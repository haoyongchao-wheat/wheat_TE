# wheat_TE
a pipeline for wheat genome TE annotation
# TE Statistics 使用示例

本文档提供了TE Statistics工具的详细使用示例，涵盖各种应用场景。

## 目录
1. [基础使用示例](#基础使用示例)
2. [小麦基因组分析示例](#小麦基因组分析示例)
3. [高级分析示例](#高级分析示例)
4. [可视化定制示例](#可视化定制示例)
5. [批量处理示例](#批量处理示例)
6. [性能优化示例](#性能优化示例)

## 基础使用示例

### 示例1: 最简单的分析
```bash
# 对GFF3文件进行基础分析
python te_statistics.py -i annotation.gff -o basic_results/
```

**输出内容:**
- `te_statistics_detailed.csv`: 详细统计表
- `plots/`: 包含5种基本图表
- `execution_summary.txt`: 执行摘要

### 示例2: 指定基因组大小
```bash
# 为计算基因组占比提供基因组大小参数
python te_statistics.py \
    -i annotation.gff \
    -o results_with_genome/ \
    -g 16000000
```

**输出特点:**
- 包含基因组占比百分比
- 每Mb转座子密度分析

### 示例3: 多格式输出
```bash
# 生成CSV和Excel格式输出
python te_statistics.py \
    -i annotation.gff \
    -o multi_format_results/ \
    --format csv excel
```

**输出文件:**
- `te_statistics_detailed.csv`
- `te_statistics_detailed.tsv`
- `te_statistics_report.xlsx` (包含多个工作表)

## 小麦基因组分析示例

### 示例4: 小麦基因组完整分析
```bash
# 使用小麦专用配置文件
python te_statistics.py \
    -i wheat_annotation.gff3 \
    -o wheat_analysis/ \
    --config te_stats/config/wheat_config.yaml \
    --genome-size 16000000000 \
    --format excel csv \
    --plot-format png svg pdf \
    --exclude-unspecified \
    --verbose
```

**配置说明:**
- `genome-size`: 16Gb (六倍体小麦)
- `exclude-unspecified`: 排除未分类TE
- `wheat_config.yaml`: 预定义的小麦TE分类规则

### 示例5: 小麦亚基因组分析
```bash
# 针对小麦特定染色体的分析
python te_statistics.py \
    -i wheat_chr1A.gff3 \
    -o wheat_chr1A_analysis/ \
    --genome-size 2700000000 \
    --min-length 200 \
    --plot-type chromosome dashboard \
    --verbose
```

## 高级分析示例

### 示例6: 自定义过滤条件
```bash
# 严格的过滤条件分析
python te_statistics.py \
    -i large_genome.gff3 \
    -o strict_analysis/ \
    --min-length 500 \
    --exclude-unspecified \
    --genome-size 3000000000 \
    --format excel \
    --verbose
```

**过滤参数:**
- `--min-length 500`: 只分析≥500bp的TE
- `--exclude-unspecified`: 排除未分类TE

### 示例7: 仅统计分析（不生成图表）
```bash
# 快速统计分析，跳过可视化
python te_statistics.py \
    -i dataset.gff3 \
    -o stats_only/ \
    --format csv excel \
    --verbose
```

### 示例8: 仅生成图表
```bash
# 基于已有数据仅生成图表
python te_statistics.py \
    -i dataset.gff3 \
    -o plots_only/ \
    --plot-only \
    --plot-type all \
    --plot-format svg pdf \
    --style ggplot
```

## 可视化定制示例

### 示例9: 特定图表类型
```bash
# 只生成饼图和柱状图
python te_statistics.py \
    -i dataset.gff3 \
    -o specific_plots/ \
    --plot-type pie bar \
    --plot-format png svg \
    --style seaborn-v0_8
```

**生成图表:**
- `composition_class1.png/svg`
- `composition_class2.png/svg`
- `comparison_class2.png/svg`

### 示例10: 高质量图表输出
```bash
# 生成出版级质量图表
python te_statistics.py \
    -i publication_data.gff3 \
    -o publication_figures/ \
    --plot-type all \
    --plot-format png svg pdf \
    --style seaborn-v0_8 \
    --verbose
```

**图表特点:**
- 300 DPI分辨率
- Arial字体
- 科学期刊级别配色
- 多种格式输出

### 示例11: 自定义样式分析
```bash
# 使用不同可视化风格
python te_statistics.py \
    -i dataset.gff3 \
    -o custom_style/ \
    --plot-type dashboard \
    --plot-format png \
    --style classic
```

**可用样式:**
- `seaborn-v0_8` (默认)
- `ggplot`
- `classic`
- `bmh`
- `dark_background`

## 批量处理示例

### 示例12: 脚本化批量分析
```bash
#!/bin/bash
# 批量分析多个文件

# 创建结果目录
mkdir -p batch_results

# 文件列表
files=("chr1A.gff" "chr1B.gff" "chr1D.gff" "chr2A.gff" "chr2B.gff" "chr2D.gff")

# 循环处理每个文件
for file in "${files[@]}"; do
    echo "Processing $file..."
    python te_statistics.py \
        -i "$file" \
        -o "batch_results/${file%.*}_analysis/" \
        --genome-size 2700000000 \
        --format csv \
        --plot-type dashboard \
        --quiet
done

echo "Batch processing completed!"
```

### 示例13: Python批量处理脚本
```python
#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

def batch_analyze_tes(input_files, output_base_dir, config_file=None):
    """
    批量分析多个GFF文件
    """
    output_dir = Path(output_base_dir)
    output_dir.mkdir(exist_ok=True)

    for gff_file in input_files:
        file_name = Path(gff_file).stem
        result_dir = output_dir / f"{file_name}_analysis"

        print(f"Processing {gff_file}...")

        cmd = [
            "python", "te_statistics.py",
            "-i", gff_file,
            "-o", str(result_dir),
            "--genome-size", "16000000",
            "--format", "csv",
            "--plot-type", "dashboard"
        ]

        if config_file:
            cmd.extend(["--config", config_file])

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✓ Completed: {file_name}")
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed: {file_name} - {e}")

# 使用示例
if __name__ == "__main__":
    gff_files = ["dataset1.gff", "dataset2.gff", "dataset3.gff"]
    batch_analyze_tes(gff_files, "batch_results")
```

## 性能优化示例

### 示例14: 大文件优化处理
```bash
# 针对大文件的优化配置
python te_statistics.py \
    -i large_genome.gff3 \
    -o large_genome_analysis/ \
    --genome-size 5000000000 \
    --min-length 100 \
    --format csv \
    --plot-type composition dashboard \
    --verbose
```

**优化策略:**
- 增加最小长度过滤减少数据量
- 选择性生成图表类型
- 使用CSV格式减少内存占用

### 示例15: 内存限制处理
```bash
# 内存受限环境下的处理
python te_statistics.py \
    -i huge_dataset.gff3 \
    -o memory_efficient_analysis/ \
    --min-length 200 \
    --exclude-unspecified \
    --format csv \
    --plot-only \
    --plot-type pie \
    --verbose
```

## 特殊应用示例

### 示例16: 比较分析
```bash
# 两个不同品系的比较分析
python te_statistics.py \
    -i cultivar_A.gff3 \
    -o cultivar_A/ \
    --genome-size 16000000000 \
    --format excel

python te_statistics.py \
    -i cultivar_B.gff3 \
    -o cultivar_B/ \
    --genome-size 16000000000 \
    --format excel
```

### 示例17: 时间序列分析
```bash
# 不同发育阶段的TE活性分析
stages=("seedling" "vegetative" "flowering" "mature")

for stage in "${stages[@]}"; do
    python te_statistics.py \
        -i "${stage}_annotation.gff3" \
        -o "time_series/${stage}/" \
        --genome-size 16000000000 \
        --plot-type composition length \
        --format csv
done
```

## 配置文件示例

### 示例18: 自定义配置文件
```yaml
# custom_analysis_config.yaml
classification:
  level1:
    "Class I": ["LTR", "LINE", "SINE", "DIRS"]
    "Class II": ["TIR", "Helitron", "Maverick", "DNA"]
    "Other": ["Unspecified", "Unknown"]

analysis:
  genome_size: 3000000000
  merge_overlaps: true
  min_length: 100
  exclude_unspecified: true

output:
  formats: ["csv", "excel"]
  create_plots: true
  plot_types: ["pie", "dashboard"]

visualization:
  style: "ggplot"
  dpi: 600
  figure_formats: ["png", "svg"]

colors:
  level1:
    "Class I": "#2E8B57"
    "Class II": "#DC143C"
    "Other": "#4682B4"
```

使用自定义配置:
```bash
python te_statistics.py \
    -i data.gff3 \
    -o custom_results/ \
    --config custom_analysis_config.yaml
```

## 结果解读示例

### 示例19: 结果文件分析
```python
import pandas as pd
import matplotlib.pyplot as plt

# 读取详细统计结果
df = pd.read_csv('results/te_statistics_detailed.csv')

# 按Level 1分组统计
level1_stats = df.groupby('level1').agg({
    'element_count': 'sum',
    'merged_length': 'sum',
    'genome_percentage': 'sum'
}).sort_values('merged_length', ascending=False)

print("Level 1 Summary:")
print(level1_stats)

# 可视化基因组占比
plt.figure(figsize=(10, 6))
level1_stats['genome_percentage'].plot(kind='pie', autopct='%1.1f%%')
plt.title('TE Genome Composition by Level 1')
plt.ylabel('')
plt.savefig('custom_genome_composition.png', dpi=300, bbox_inches='tight')
plt.show()
```

## 故障排除示例

### 示例20: 调试模式运行
```bash
# 详细调试信息
python te_statistics.py \
    -i problematic_file.gff3 \
    -o debug_results/ \
    --verbose 2>&1 | tee debug.log

# 检查日志文件
tail -f te_statistics.log
```

### 示例21: 分步验证
```bash
# 第1步：仅验证GFF解析
python -c "
from te_stats.src.gff_parser import GFFParser
parser = GFFParser('test.gff')
records = parser.parse_gff_file()
print(f'Parsed {len(records)} records')
"

# 第2步：验证分类器
python -c "
from te_stats.src.gff_parser import GFFParser
from te_stats.src.classifier import TEClassifier
parser = GFFParser('test.gff')
records = parser.parse_gff_file()
classifier = TEClassifier()
for i, record in enumerate(records[:10]):
    print(f'{record.te_class} -> {classifier.classify_element(record.te_class)}')
"
```

## 最佳实践建议

1. **数据准备**: 确保GFF3文件格式正确，包含必要的Class字段
2. **内存管理**: 大文件建议增加最小长度过滤阈值
3. **输出选择**: 根据需要选择合适的输出格式和图表类型
4. **配置管理**: 为不同物种创建专用的配置文件
5. **结果验证**: 始终检查执行摘要和日志文件
6. **批量处理**: 使用脚本进行批量分析的标准化

通过这些示例，您可以根据具体需求灵活使用TE Statistics工具进行转座子分析。
