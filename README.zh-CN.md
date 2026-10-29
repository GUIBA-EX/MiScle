# MiScle

[English](README.md)

**MiScle, a set of PCR primer candidates for metabarcoding environmental DNA from scleractinian corals**

## 摘要

MiScle比较石珊瑚线粒体16S的S、M两个同源窗口。S保留普通与D引物；M采用一条共用F和Refertina／Vacatina两个R，普通方案为MR＋MV，简并方案为MRD＋MVD。R/V表示优先补偿的亚目，D表示简并；两个亚目均为目标，允许交叉覆盖。每种M方案只需合成三种寡核苷酸，共用F只加入一次。

在M窗口304个可评价物种中，普通组合预测覆盖282个（92.76%）、完全匹配78个（25.66%）；D组合覆盖284个（93.42%）、完全匹配98个（32.24%）。原M和MD均覆盖124个，分别完全匹配1和81个。新组合的严格WC/GT背景并集均为103/111条，原方案为0/0。它们是覆盖优先的探索性子引物，增加背景命中和Vacatina完全匹配下降是实质取舍。

共用F条件下，M-D组合在Vacatina中完全匹配12/154（7.79%），低于原MD的75/154（48.70%）。S与M-D合并后完全匹配160/305，也低于原SD＋MD的206/305。因此本轮提高了M窗口的覆盖范围，尚未同时解决两个亚目的完全匹配。后续组织DNA和混合样本实验应将原M/MD作为对照。

## 引物序列与名称

所有靶向序列按合成方向5′→3′书写，不含接头。Y=C/T，R=A/G。S/M表示同源窗口，实际产物长度随谱系而异。D每条最多8倍；本次M-D三种寡核苷酸均为2倍。

| 名称 | F，5′→3′ | R，5′→3′ | F/R简并度 |
|---|---|---|---:|
| MiScle-S | `TTTGATAAGTGGGACAGTTTGG` | `TTAACGGGTCACTTAAACCC` | 1 / 1 |
| MiScle-SD | `YTKGATAAGTGRGACAGTTTGG` | `YYAACGGGTCACTTWAACCC` | 8 / 8 |
| MiScle-MR | `AAGTCCCCATGGAACTTTAC` | `CTTAAGCCCGCCTTCGTTAC` | 1 / 1 |
| MiScle-MV | `AAGTCCCCATGGAACTTTAC` | `TCATAGTTCGTTTTCGTTAC` | 1 / 1 |
| MiScle-MRD | `AAGTCCCYATGGAACTTTAC` | `CTTAAGYCCGCCTTCGTTAC` | 2 / 2 |
| MiScle-MVD | `AAGTCCCYATGGAACTTTAC` | `CCATAGTTCRTTTTCGTTAC` | 2 / 2 |

普通MR与MV的F完全相同，MRD与MVD的F完全相同。四个M配对的参考端点均为F1561–1580、R1741–1760；普通与D的对应坐标完全一致。F/R的两个数分别表示每条寡核苷酸的简并度，不是整管混合物的总简并度。

[序列表](primer_pairs.csv) · [FASTA](primers.fasta) · [全部参数](selected.json) · [原方案对照](results/baseline/primer_pairs.csv)

## 材料与方法

### 参考数据与分类

Hisata等（2025）研究的是12S系统；本项目提取其表1中有注释的线粒体16S，再补充Scleractinia公开记录，包括非造礁和深海类群。110个论文accession中65条含可用16S，45条没有可提取的16S，不计作扩增失败。目标检索返回1974条记录；合并审核保留879条合格16S，以来源taxid和具名双名记录去重，每物种一条代表，优先完整基因、长度和低模糊度，不按引物匹配挑选代表。

冻结面板含307个物种标签、131个属、29个来源科级标签。亚目来源lineage为Refertina 148、Vacatina 157、未定2；框架参考WoRMS和Vaga等（2025）。Guynia和Stenocyathus仍按来源未定处理。来源科名保留旧归属，科级标签不等同于29个确认的现行有效科；分析同时报告属和物种。仅提取方向正确的16S／线粒体大亚基rRNA注释，排除短于150 nt或模糊碱基超过5%的记录。MAFFT用于多序列比对（Katoh & Standley, 2013）。

坐标为Acropora tenuis NC_003522.1提取后、去缺口16S基因的1-based位置，不是整个线粒体坐标。S为F1696–1717/R1826–1845，参考跨度150 bp；M为F1561–1580/R1741–1760，参考跨度200 bp。覆盖参考中的实际含引物产物为S 147–209 bp、MR/MRD 186–298 bp、MV/MVD 99–186 bp。

![参考结合位置](figures/reference_windows.svg)

### 同坐标亚目子引物搜索

原设计参考MiFish的保守锚点、谱系补偿和成对评价思路（Miya et al., 2015）。本轮固定M四个端点，使用独立复核的原始结合模板及未命中物种的同源列生成候选。分别构建共用F、Refertina R、Vacatina R候选池。普通候选使用真实ACGT模板及共识；D参考Anth-28S-eDNA的实际变异补偿逻辑（McCartin et al., 2024），迭代合并同源列变异，单条简并度≤8，不限制仅在3′半段简并。

保留40个F、13个Refertina R和32个Vacatina R候选，去重后79条寡核苷酸。全部进行双方向原始序列位点扫描，再按成对方向、间距和末端兼容性评价。比较低背景筛选与放宽背景筛选；原低背景条件没有得到充分的Refertina互补收益。当前候选的探索性筛选上限为每配对严格WC≤150、GT≤150、WC3≤1200条背景；该上限不表示背景扩增可忽略。在最佳预测覆盖至多损失2个物种的范围内，优先完全匹配，其次覆盖和较弱亚目的覆盖。搜索为有限候选的启发式比较，不保证全局最优。

### 热力学与混合物

Primer3采用单价盐50 mM、Mg²⁺1.5 mM、dNTP 0.6 mM（Untergasser et al., 2012）。比较每具体变体50 nM及共用F总量200 nM、两个R各100 nM两种情景；D三条各2倍时，后一情景每F变体100 nM、每R变体50 nM。检查全部展开序列的发卡、自二聚体和交叉二聚体，包括两个R之间的互作。要求最高结构Tm<55°C、F与任一R的双链Tm区间中点差≤5°C、最低双链Tm≥45°C。两种M混合物最高预测结构Tm均为44.46°C。末端3 nt至少有确定的G/C，排除可能形成五连同聚物的候选。加尾全长尚未评价。

### 覆盖、背景和独立复核

覆盖允许每条≤2个错配，末端3 nt按IUPAC兼容；完全匹配要求两端均零错配。原始序列扫描检查双方向，目标产物70–500 bp且内部≥20 bp；背景内部间距0–1000 bp。WC为Watson–Crick模型，GT允许计算性G–T配对；WC3、GT2、WC4为宽松敏感性检查，均不等同于PCR检出概率。

背景为9145条冻结记录：1938条分层刺胞动物参考和7207条棘皮动物参考。近缘背景有查询上限，未穷尽所有海洋生物，也未系统覆盖细菌与共生藻。M有304个可评价物种、3个未知；S和M至少一个可评价的并集为305。缺少位点信息不当作失败。

六个正式配对与八个额外混合配对（每种M方案的F/F、R/R、V/V、R/V）均完成独立C++穷举和Python位并行扫描，目标及严格背景覆盖／完全匹配共84项集合比较，差异为0。八个额外配对在目标及严格WC/GT背景中均无命中。本轮同管指标是序列兼容性并集及混合物互作的计算检查，不能证明相同PCR条件下没有竞争或扩增偏倚。跨S、M组合按分管检测后的分类结果并集解释。

## 预测结果

| 方案 | 可评价 | 覆盖 | 完全匹配 | 严格WC/GT | 宽松WC3/GT2/WC4 |
|---|---:|---:|---:|---:|---:|
| MiScle-S | 304 | 153 (50.33%) | 64 (21.05%) | 0 / 0 | 0 / 0 / 0 |
| MiScle-SD | 304 | 155 (50.99%) | 125 (41.12%) | 0 / 0 | 6 / 0 / 211 |
| MiScle-MR | 304 | 160 (52.63%) | 78 (25.66%) | 103 / 111 | 542 / 111 / 1152 |
| MiScle-MV | 304 | 122 (40.13%) | 0 (0.00%) | 0 / 0 | 0 / 0 / 0 |
| MiScle-MRD | 304 | 160 (52.63%) | 93 (30.59%) | 103 / 111 | 645 / 111 / 1168 |
| MiScle-MVD | 304 | 124 (40.79%) | 5 (1.64%) | 0 / 0 | 0 / 0 / 4 |

### M组合与原方案比较

| 方案 | 可评价 | 覆盖 | 完全匹配 | 严格WC/GT背景 |
|---|---:|---:|---:|---:|
| 原M | 304 | 124（40.79%） | 1（0.33%） | 0 / 0 |
| MR＋MV | 304 | 282（92.76%） | 78（25.66%） | 103 / 111 |
| 原MD | 304 | 124（40.79%） | 81（26.64%） | 0 / 0 |
| MRD＋MVD | 304 | 284（93.42%） | 98（32.24%） | 103 / 111 |

M-D组合分别覆盖Refertina 145/148（97.97%）、Vacatina 137/154（88.96%），完全匹配85/148（57.43%）和12/154（7.79%）；两个未定物种均覆盖，其中1个完全匹配。M-D仍漏检20个可评价物种。跨窗口SD＋MRD＋MVD覆盖285/305（93.44%）、完全匹配160/305（52.46%）；原SD＋MD对应为279/305（91.48%）和206/305（67.54%）。

103条WC背景全部来自Actiniaria（70）和Corallimorpharia（33）；GT为72和39条。按全部背景分母为1.13%和1.21%，按1938条刺胞动物参考分母则为5.31%和5.73%。Refertina子引物承担了这些新增命中，不能概括为仍保持原来的零背景。MRD的WC3/WC4为645/1168条，需在近缘非目标混样中检验。

[逐物种](by_species.csv) · [逐亚目](by_suborder.csv) · [亚目组合](combination_by_suborder.csv) · [组合比较](two_pair_comparison.csv) · [背景命中清单](background_hits.csv) · [M-D漏检](missed_species.csv) · [原方案结果](results/baseline/comparison.csv)

## 碱基匹配与条码分辨率

图按亚目、来源科分组，两端均按5′→3′显示，R模板已取反向互补。E为完全匹配，C为允许错配下覆盖，–为无预测命中，?为不可评价。蓝色为兼容、橙色为错配、灰色为缺口／未知，方框标出末端3 nt。覆盖行采用独立复核的原始结合模板，其他行采用参考同源列。R/V名称说明设计侧重，不表示只扩增对应亚目。

[碱基数据](base_matching.csv) · [实际结合模板与产物](verified_templates.csv) · [条码分辨率](barcode_resolution.csv) · [共享条码](shared_barcodes.csv)

MR/MRD中48/160个覆盖物种的内部条码未与其他覆盖物种共享；MV为24/122，MVD为26/124。共享条码和跨属共享仍存在。这是当前参考库的观察，不是未知样品鉴定成功率；结合位点改善不会自动提高内部条码的物种分辨率。

![MiScle-S](figures/base_matching_MiScle-S.svg)

![MiScle-SD](figures/base_matching_MiScle-SD.svg)

![MiScle-MR](figures/base_matching_MiScle-MR.svg)

![MiScle-MV](figures/base_matching_MiScle-MV.svg)

![MiScle-MRD](figures/base_matching_MiScle-MRD.svg)

![MiScle-MVD](figures/base_matching_MiScle-MVD.svg)

## 复现与文件

- `python3 scripts/check_release.py`：核对名称、序列、坐标、共用F、简并度、图和结果。
- `python3 scripts/verify.py`：使用冻结输入独立重扫，Python标准库即可，无需联网。
- [复核摘要](verification.json) · [84项比较](independent_checks.csv) · [搜索及取舍](results/suborder_search/README.md) · [参数与输入来源](manifest.json) · [accession审核](data/paper_accession_audit.csv)。
- 新的候选生成和配对筛选脚本见上述搜索说明；其他早期设计脚本保留工作区路径，作为过程记录。最终文件校验和见 `SHA256SUMS`。

当前优先实验是原M/MD与新R/V方案的组织DNA对照，以及含海葵、Corallimorpharia近缘非目标的平衡和不平衡混样。若目标是进一步提高Vacatina完全匹配，下一轮应比较放开共用F约束的亚目完整配对，而不是仅继续增加R简并。

## 参考文献

1. Hisata, K., Nagata, T., Kanai, M., et al. (2025). An eDNA metabarcoding system for detecting scleractinian corals to the generic level along the Japanese coast. *Galaxea, Journal of Coral Reef Studies*, 27, 13–29. https://doi.org/10.3755/galaxea.G27D-5
2. Miya, M., Sato, Y., Fukunaga, T., et al. (2015). MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species. *Royal Society Open Science*, 2, 150088. https://doi.org/10.1098/rsos.150088
3. McCartin, L., Saso, E., Vohsen, S. A., et al. (2024). Nuclear eDNA metabarcoding primers for anthozoan coral biodiversity assessment. *PeerJ*, 12, e18607. https://doi.org/10.7717/peerj.18607
4. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30, 772–780. https://doi.org/10.1093/molbev/mst010
5. Untergasser, A., Cutcutache, I., Koressaar, T., et al. (2012). Primer3—new capabilities and interfaces. *Nucleic Acids Research*, 40, e115. https://doi.org/10.1093/nar/gks596
6. Vaga, C. F., Quattrini, A. M., Galvão de Lossio e Seiblitz, I., et al. (2025). A global coral phylogeny reveals resilience and vulnerability through deep time. *Nature*, 648, 377–382. https://doi.org/10.1038/s41586-025-09615-6
7. World List of Scleractinia. [Taxonomic classification](https://www.marinespecies.org/scleractinia/aphia.php?p=stats). Accessed 2026-10-30.

## 著作权

Copyright © 2026 Fei Xia (GUIBA-EX). 保留所有权利。未经书面许可禁止转载和再分发受保护的原创材料。公共数据库记录保留原有权利状态，详见[著作权与使用声明](COPYRIGHT.md)。
