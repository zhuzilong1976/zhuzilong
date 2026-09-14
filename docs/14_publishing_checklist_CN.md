# 发布清单：GitHub 仓库 + Zenodo DOI

更新日期：2026-09-14 ｜ 对应 `docs/09_methods_and_limitations_draft_EN.md` 的
`Data and code availability` 与根目录 `CITATION.cff`

---

## 0. 四个值的落实情况（2026-09-14 已填）

| 信息 | 值 | 状态 |
|---|---|---|
| 作者 + ORCID | 朱子龙 / **Zhu Zilong**，ORCID `0000-0002-6955-0903` | ✅ 已写入 `CITATION.cff` |
| GitHub 仓库地址 | `https://github.com/zhuzilong1976/zhuzilong` | ✅ 2026-09-14 已确认存在、**public**、`default_branch = main`；当前仍是空仓库，本地 5 个提交待推送 |
| 论文题目 | `Containment, not biology: the differentially regulated gene set of single-cell virtual knockout is the knockout's own out-neighbourhood` | ✅ 已填入 `preferred-citation`（暂定） |
| 目标期刊 | **`Briefings in Bioinformatics`**（2026-09-14 确定） | ✅ 已写入 `preferred-citation` |
| 文章类型 | **Original Article**（2026-09-14 确定） | ✅ 已写入 cover letter 第 1 段 |
| 作者单位 | **Department of Neurology, Tianjin Huanhu Hospital, Tianjin, China**（天津市环湖医院神经内科） | ✅ 已写入正文标题页与 cover letter 署名 |
| git 提交身份 | `zhuzilong` / `zhuzilong1976@gmail.com` | ✅ 已配置为仓库级 `user.name`/`user.email` 并完成首次提交 |

只剩 Zenodo DOI 一项（第 3 节）。

---

## 1. 本地仓库状态

仓库已在 `vko-project/` 下初始化为 git 仓库，并已完成首次提交：

```text
37b0cc3  Initial release: code, figures and tables for the virtual-knockout containment study
         109 files, 0 files > 1 MB, working tree clean
```

尚未 `push`（需要你的 GitHub 凭据）。

`.gitignore` 已排除大文件与派生数据：`work/`、`outputs/`、`data/*`、`*.rds`、
`*.h5ad`、`*.pkl`、`*.mtx`、`*.zip`。因此入 git 的只有代码、文档、图与表
（约 1 MB）；完整距离矩阵等大文件走 Zenodo（第 3 节）。

## 2. 建立 GitHub 仓库并推送

```bash
cd C:\Users\朱子龙\Desktop\codex\vko-project

git branch -M main

# 先在 GitHub 网页新建一个空仓库（不要勾选 README/.gitignore），然后：
git remote add origin https://github.com/zhuzilong1976/zhuzilong.git
git push -u origin main
```

推送前用 `git status` 确认没有 `work/`、`outputs/` 下的文件。

## 3. Zenodo 归档拿 DOI（推荐 GitHub 集成）

1. 登录 Zenodo → Settings → GitHub → 打开该仓库的开关；
2. 回 GitHub 创建一个 Release（例如 tag `v1.0.0`），Zenodo 自动归档并分配 DOI；
3. 在该页面可拿到 DOI 徽标与 concept DOI。

完整敲除距离矩阵、20 个稳定性重复对象、CellOracle 结果不随 git 走，两种处理方式：

- **方式 A（推荐）**：这些数据可由 `code/` 重新生成，`data/README.md` 已写明输入数据来源；
  `Data and code availability` 按"代码在 GitHub、原始数据在 GEO"表述（当前 `docs/09` 即如此）；
- **方式 B**：想连派生数据一起归档，先打包再上传 Zenodo：

  ```powershell
  powershell -File code\00_setup\make_transfer_bundle.ps1 -Dest D:\transfer -Full
  Compress-Archive -Path D:\transfer\vko-project -DestinationPath D:\transfer\vko-project-full.zip
  ```

  然后把该 zip 作为 Zenodo 的一个版本上传（GitHub 单文件上限 100 MB，不要放进 git）。

## 4. 拿到 DOI 之后回填三处

1. `docs/09_methods_and_limitations_draft_EN.md` 的 `Data and code availability`：
   把 `**[待补: DOI]**` 换成真实 DOI；
2. `CITATION.cff` 增加 `doi: "10.5281/zenodo.XXXXXXX"`，并补上作者、`repository-code`、
   `preferred-citation`；
3. `README.md` 的 Citation 段可加 DOI 徽标。

## 7. 目标期刊：Briefings in Bioinformatics（2026-09-14 确定）

### 7.1 已核实的指标（2026 年 6 月发布的 JCR 2025 数据）

| 期刊 | IF 2025 | IF 2024 | IF 2023 | 中科院大类分区 | 大类学科 |
|---|---:|---:|---:|---|---|
| **Briefings in Bioinformatics** ← 目标 | **7.3** | 7.7 | 6.8 | **2区** | 生物学（小类：生化研究方法） |
| Bioinformatics (OUP) | 5.5 | 5.4 | 4.4 | 2区 | 生物学（生化研究方法） |
| PLOS Computational Biology | 3.7 | 3.6 | 3.8 | 2区 | 生物学（生化研究方法） |
| Genome Biology | 9.2 | 9.4 | 10.1 | 1区 | 生物学（生物工程与应用微生物） |
| Nature Methods | 28.3 | 32.1 | 36.1 | 1区 | 生物学（生化研究方法） |
| NAR Genomics and Bioinformatics | 3.0 | – | – | 3区 | 生物学（遗传学） |
| GigaScience | 5.0 | – | – | 3区 | 生物学（综合性期刊） |
| Cell Reports Methods | 5.8 | 4.5 | – | 2区 | 生物学（生化研究方法） |

来源：BioxBio（历年 IF）与 LetPub 期刊库（IF、中科院分区）。两者对 BiB 的 2025 年 IF 一致（7.3）。
JCR 类别排名（Q1/Q2）请以 Clarivate JCR 官方页面为准；OUP 官网目前被 Cloudflare 拦截，
因此**文章类型、字数上限与 APC 需在投稿系统内确认**。

### 7.2 BiB 的其他参考信息（LetPub 网友数据，供风险评估）

- 自引率 6.8%；CiteScore 13.60；SJR 2.264；h-index 90
- 审稿速度：平均约 6 个月；录用难度：较难；Gold OA：No（混合期刊，需确认是否收取版面费）

### 7.3 投稿前必须确认的四件事

1. **文章类型**：已定为 **Original Article**（机制刻画 + 原始结果 + 可复现性评估 + 实践建议）。
   投稿时请对齐投稿系统里的实际标签（可能写作 Research Article / Original Research），若不同只需改
   cover letter 第 1 段那一处；
2. **APC 与版权**：混合期刊，选择 CC BY 开放获取才产生费用，需确认；
3. **字数与图表限制**：Results 目前含 Table S10/S11 与 Figure S2，正文篇幅需按新限制裁剪；
4. **Cover letter**：需在开头一段讲清"这是对已发表工具的方法学评估，并给出可操作的报告规范建议"。

## 5. 归档前自查

- [x] `code/` 与实测一致：本轮只**新增** `06_figures/54_figure_stability.R` 与
      `07_supplementary/67_stability_checkpoints_C3C4.R`，既有脚本一行未改；
- [x] 新增脚本可复现：`67_…` 归档后重跑，两个 checkpoint CSV 与归档前**逐字节一致**；
- [x] `results/figures/` 含 Figure 1–3、FigS1、**FigS2**；`results/tables/` 含 S1–**S11**；
- [x] README 的计数、图表范围、复现步骤、已知缺口已与实测一致（38 网络 / 约 20,800 次敲除）；
- [x] LICENSE：代码 MIT、`results/` 下图表 CC BY 4.0、第三方数据来源已列明；
- [ ] `CITATION.cff` 的 TBD（作者、repo、preferred-citation）——需要第 0 节的信息；
- [ ] Zenodo DOI——需要第 3 节操作；
- [ ] 若期刊要求"代码可运行性"证明，可把 `docs/12_runtime_benchmarks_CN.md` 与
      `vko-run/COMPLETION_REPORT_CN.md` 第 6 节的本机实测合并成一段补充说明。

发布前顺手修掉的三处问题（本轮已改）：

- `logs/`（本机稳定性运行日志、含绝对路径）加入 `.gitignore`，不再入库；
- `docs/08_cross_dataset_replication_CN.md` 里指向原机器绝对路径的图片链接改为相对路径
  `../results/figures/Fig_replication.png`；
- `env/python_requirements.txt` 顶部加了说明：该文件里多数条目是原机器的本地 wheel 路径，
  不能直接 `pip install -r`，版本号在 wheel 文件名里。

## 6. 本轮（2026-09-14）已完成

| 项目 | 位置 |
|---|---|
| C3/C4 检查点脚本归档 | `code/07_supplementary/67_stability_checkpoints_C3C4.R` |
| 稳定性四联图 | `code/06_figures/54_figure_stability.R` → `results/figures/FigS2_stability.{png,pdf}` |
| 稳定性逐重复表与判定表 | `results/tables/Table_S10_stability_replicates.csv`、`Table_S11_stability_summary.csv` |
| 深度匹配结果 | `outputs/results/vko_ms_depthmatched/`（27 个文件，含 `panel_depthmatched_vs_main.csv`） |
| 稿件 Methods/Results/Limitations 改写 | `docs/09_methods_and_limitations_draft_EN.md`（原稿备份在 `vko-run/backup/`） |
| 运行记录与运维脚本 | `vko-run/`（`COMPLETION_REPORT_CN.md`、`README_operations_CN.md`、日志） |
