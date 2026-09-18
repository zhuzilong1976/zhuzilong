# Cover letter — Briefings in Bioinformatics

> Targets the manuscript in `docs/15_manuscript_full_draft_EN.md`. Drafted to the
> `cover-letter-drafter` specification: five paragraphs, 300–450 words, no claims about
> acceptance likelihood, no invented affiliations or editorial preferences.
> Placeholders are marked `[TO BE COMPLETED]` / `[Author to confirm]`.

---

Dear Editors of *Briefings in Bioinformatics*,

**Submission.** We submit for consideration the manuscript "Containment, not biology: the differentially regulated gene set of single-cell virtual knockout is the knockout's own out-neighbourhood", as a Problem Solving Protocol. The journal reaches the readers who run single-cell gene regulatory network analyses, and the manuscript answers a question they meet whenever they interpret a virtual-knockout gene list: what does that list contain, and what must be reported with it? We diagnose the problem for a tool in active use (`scTenifoldKnk`) with results across 3,195 knockouts and 10 published networks, and set out the reporting protocol that follows: state the null model, report the outdegree of every knockout, do not rank genes by list size, and pair every condition comparison with a same-condition null.

**What is new.** Benchmarking has shown that in-silico perturbation results are determined by method choice rather than biology (Wu et al., 2026), and confounder catalogues are appearing (Qiu and Zhao, 2026); the mechanism behind an individual tool had not been identified. The differentially regulated (DR) set is not a downstream readout but a structural property of the inferred network: it is contained in the knocked-out gene plus its direct out-neighbours, and its size is set by the knockout's outdegree through dilution of each target's incoming weight.

**Approach and principal results.** From 3,195 knockouts in four human MS-lesion microglia networks, replicated in 10 published networks (2,591–9,970 genes; human and mouse) and in the developers' own *Trem2* result, containment held in 99.87% of knockouts and in 100% of non-degenerate knockouts. Two consequences matter in practice. The default χ² null returns exactly one significant gene — the knocked-out gene itself — in 96% of knockouts, whereas the empirical null returns outdegree + 1, so the null model must be stated. Biologically irrelevant knockouts can also dominate: in the MS network a negative-control gene with outdegree 0 produced 188 significant genes (adjusted *p* down to 1.3 × 10⁻²¹⁶), while the same gene returned none in the other condition. The property is tool-specific: with CellOracle, 44.5–77.9% of affected genes lie outside the direct target set.

**Relevance to your readers.** Four requirements follow: report the outdegree of every knockout; do not rank genes by the number of significant targets; pair every condition comparison with a same-condition null; and state the null model used. The null comparison is decisive — the MS-versus-control rate (7.25% at *p* < 0.05) was indistinguishable from a split-half null of the same condition (5.65%; Fisher *p* = 0.221), and the null's top genes were equally plausible, so a convincing list is not evidence. Code, network objects and all knockout distance matrices are openly available (https://github.com/zhuzilong1976/zhuzilong; archived at Zenodo, DOI 10.5281/zenodo.22752977).

**Declarations.** The manuscript is original, has not been published and is not under consideration elsewhere; all authors have approved submission. The authors declare no competing interests. Thank you for considering this work.

Sincerely,

Zilong Zhu
ORCID: 0000-0002-6955-0903
Department of Neurology, Tianjin Huanhu Hospital, Tianjin, China
zhuzilong1976@gmail.com

---

## Optional block: suggested reviewers

Insert into the submission form if the journal requests reviewers. Names, institutions, ORCIDs and e-mail addresses below were taken from the corresponding-author line of an open-access full text (PMC) or from the publisher's article page, and were verified on 2026-09-17; no address was guessed.

| # | Candidate (verified) | Basis for expertise | Rationale |
|---|---|---|---|
| 1 | **Prof. T. M. Murali**, Virginia Tech, USA ｜ **murali@cs.vt.edu** ｜ ORCID 0000-0003-3688-4672 | Last author of "Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data", *Nature Methods* 2020 (doi:10.1038/s41592-019-0690-6; ~870 citations) | Senior authority on how single-cell GRN inference should be evaluated; can judge whether the containment rule generalises beyond the tested tools |
| 2 | **Prof. Zhao-Peng Liu**, Shandong University, China ｜ **zpliu@sdu.edu.cn** | Author of single-cell GRN-inference methods and of the corresponding evaluations in *Briefings in Bioinformatics* (GeneLink 2025, doi:10.1093/bib/bbaf359; LogicSR 2025, doi:10.1093/bib/bbaf621) | Independent of both tools compared here; works directly on single-cell GRN inference and its evaluation |
| 3 | **Prof. Samantha A. Morris**, Washington University in St. Louis, USA ｜ **s.morris@wustl.edu** | Senior and corresponding author of CellOracle (*Nature* 2023, doi:10.1038/s41586-022-05688-9) | Expert on the alternative, multi-hop perturbation design used as the contrast; **note: her laboratory developed a tool compared in this work** |

Further alternates if the editor wants additional names:

| # | Candidate (verified) | Basis for expertise | Rationale |
|---|---|---|---|
| 4 | **Prof. Jianlin Cheng**, University of Missouri, USA ｜ **chengji@missouri.edu** | Senior author of "Machine learning methods for gene regulatory network inference" (*Briefings in Bioinformatics* 2025) | Independent; machine-learning perspective on GRN inference |
| 5 | **Dr. Wenjie Ge / Dr. Shixiang Wu**, Wuxi Taihu Hospital, China | Systematic evaluation of eight in-silico perturbation methods across four datasets (2026 preprint, doi:10.64898/2026.08.11.744106) | Closest match to the topic; e-mail is shown on the preprint's correspondence line (not machine-readable from here) |
| 6 | **Dr. Ru Qiu**, Sun Yat-sen University, China ｜ ORCID 0009-0007-0021-7241 | Author of a confound-diagnostic toolkit for in-silico perturbation (2026 preprint, doi:10.64898/2026.08.04.732812) | Independent; works specifically on confounders of perturbation outputs |
| 7 | **Dr. Daniel Osorio**, Texas A&M University, USA | First author of scTenifoldKnk (*Patterns* 2022, doi:10.1016/j.patter.2022.100434) | Tool developer's view; include only if the editor asks for it, and declare the relationship |

E-mail addresses are listed for candidates 1–4. Candidate 5 (Ge/Wu) has no machine-readable
address: copy it from the correspondence line on the preprint page if you prefer that pair.

Reviewers who authored `scTenifoldKnk` are deliberately **not** among the three primary suggestions (candidate 5 is offered only as an alternate, for the editor's discretion), to keep the assessment independent.

---

## Internal notes (not for submission)

* Letter body: 448 words (within the 300–450 range for an Original Article submission). If the journal imposes a shorter limit, cut paragraph 3 to its first and third sentences (~330 words).
* **Article type: Original Article** (author decision, 2026-09-14). When submitting, match the exact label used in the journal's submission system; if the closest label is "Research Article" or "Original Research", change the phrase "as an Original Article" in paragraph 1 only.
* The letter does not state or imply acceptance likelihood, acceptance rates or editorial preferences.
* If the author list grows, replace "the author" with "all authors" in the Declarations paragraph and add a contributions sentence.
