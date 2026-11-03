# Brain Atlases & Label Dictionaries

日常分析常用的脑区模板（atlas）与标签字典（label dictionary）汇总，主要用于 SPM / DPABI / CAT12 / REST 等工具的 ROI 提取与结果可视化。

## 目录内容

### `AAL3/` — AAL3 图谱（ Automated Anatomical Labeling 3 ）

| 文件 | 说明 |
|---|---|
| `AAL3v1.nii`、`AAL3v1_1mm.nii.gz` | AAL3v1 分区模板（MNI 空间，2 mm / 1 mm 两种分辨率） |
| `AAL3v1.nii.txt`、`AAL3v1_1mm.nii.txt` | 每个 ROI 编号对应的英文名称 |
| `ROI_MNI_V7.nii`、`ROI_MNI_V7_1mm.nii` | ROI_MNI_V7 系列模板（含对应的 `_List.mat`、`_vol.mat`、`_Border.mat`） |
| `AAL3_ROI_124_retained.txt` | 保留 ROI 列表 |
| `AAL3_UserGuide_April2024.pdf` | 官方使用说明 |

**中英字典（自制）**

| 文件 | 说明 |
|---|---|
| `AAL3中英文字典.txt` | 主字典，TSV 格式，含 170 条对照 |
| `AAL3中英文字典.xlsx` | 同上，Excel 版本 |
| `AAL3中英文字典_筛选150ROI.xlsx` | 筛选后保留 150 个 ROI 的版本 |

主字典为 Tab 分隔，四列：

```
ROI_ID    English_Name    Chinese_Name    Atlas_Status
1         Precentral_L    左中央前回      有效标签
2         Precentral_R    右中央前回      有效标签
```

`Atlas_Status` 列用于标记该标签在图谱中是否有效，便于批量提取时过滤。

配套 MATLAB 脚本：`AAL3.m`、`batch_extract_AAL3_CBF.m`，以及 `gin_*.m` 系列（标签与簇的读取/显示工具）。

### `LPBA40SPM5/` — LPBA40 概率图谱（SPM5 版本）

LONI Probabilistic Brain Atlas，基于 40 名被试手工标注构建的概率图谱，SPM5 预处理流程。

| 子目录 | 内容 |
|---|---|
| `PDF/` | 56 个脑区的概率图，每区一个 `.nii.gz`（如 `L.hippocampus`、`R.precentral_gyrus`） |
| `PDFGM/` | 灰质受限版本的概率图 |
| `tissue/` | 组织概率图：`gm` / `wm` / `csf` |
| `maxprob/` | 最大概率标注图 |
| `avg/` | 40 名被试的平均 T1 模板 |
| 根目录 | `Atlas_Distribution_Agreement_v2.pdf`（标注一致性与分发说明）、重采样后的灰质标注图 |

## 使用提示

- 所有模板位于 **MNI 空间**，使用前确认你的数据已完成空间标准化。
- 提取 ROI 时建议用 `AAL3中英文字典.txt` 做编号到中文名的映射，`Atlas_Status` 为「有效标签」的行才参与统计。
- LPBA40 的 `PDF/` 是概率图（0–1 连续值），不是硬分割；需要硬分割请用 `maxprob/`。

## 来源与引用

- **AAL3**：Rolls, E. T., Huang, C. C., Lin, C. P., Feng, J., & Joliot, M. (2020). Automated anatomical labelling atlas 3. *NeuroImage*, 206, 116189.
- **LPBA40**：Shattuck, D. W., Mirza, M., Adisetiyo, V., et al. (2008). Construction of a 3D probabilistic atlas of human cortical structures. *NeuroImage*, 39(3), 1064-1080.（数据由 LONI / USC Laboratory of Neuro Imaging 分发）

请在使用时引用上述原始文献。中英对照字典为本人整理，按仓库 License 使用。
