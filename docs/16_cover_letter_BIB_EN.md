# Cover letter — Briefings in Bioinformatics

> Targets the manuscript in `docs/15_manuscript_full_draft_EN.md`. Drafted to the
> `cover-letter-drafter` specification: five paragraphs, 300–450 words, no claims about
> acceptance likelihood, no invented affiliations or editorial preferences.
> Placeholders are marked `[TO BE COMPLETED]` / `[Author to confirm]`.

---

Dear Editors of *Briefings in Bioinformatics*,

**Submission.** We submit for consideration the manuscript "Containment, not biology: the differentially regulated gene set of single-cell virtual knockout is the knockout's own out-neighbourhood", as an Original Article. The journal reaches the readers who run single-cell gene regulatory network analyses, and the manuscript answers the question they face when interpreting a virtual-knockout gene list: what does that list contain? It is a mechanistic assessment of a tool in active use (`scTenifoldKnk`), with original results across 3,195 knockouts and 10 published networks, condensed into four reporting requirements.

**What is new.** Benchmarking has shown that in-silico perturbation results are determined by method choice rather than biology (Wu et al., 2026), and confounder catalogues are appearing (Qiu and Zhao, 2026); the mechanism behind an individual tool had not been identified. The differentially regulated (DR) set is not a downstream readout but a structural property of the inferred network: it is contained in the knocked-out gene plus its direct out-neighbours, and its size is set by the knockout's outdegree through dilution of each target's incoming weight.

**Approach and principal results.** From 3,195 knockouts in four human MS-lesion microglia networks, replicated in 10 published networks (2,591–9,970 genes; human and mouse) and in the developers' own *Trem2* result, containment held in 99.87% of knockouts and in 100% of non-degenerate knockouts. Two consequences matter in practice. The default χ² null returns exactly one significant gene — the knocked-out gene itself — in 96% of knockouts, whereas the empirical null returns outdegree + 1, so the null model must be stated. Biologically irrelevant knockouts can also dominate: in the MS network a negative-control gene with outdegree 0 produced 188 significant genes (adjusted *p* down to 1.3 × 10⁻²¹⁶), while the same gene returned none in the other condition. The property is tool-specific: with CellOracle, 44–78% of affected genes lie outside the direct target set.

**Relevance to your readers.** Four requirements follow: report the outdegree of every knockout; do not rank genes by the number of significant targets; pair every condition comparison with a same-condition null; and state the null model used. The null comparison is decisive — the MS-versus-control rate (7.25% at *p* < 0.05) was indistinguishable from a split-half null of the same condition (5.65%; Fisher *p* = 0.221), and the null's top genes were equally plausible, so a convincing list is not evidence. Code, network objects and all knockout distance matrices are openly available (GitHub: `[TO BE COMPLETED]`; Zenodo: `[TO BE COMPLETED: DOI]`).

**Declarations.** The manuscript is original, has not been published and is not under consideration elsewhere; the author has approved submission. The author declares no competing interests. `[Author to confirm]` Thank you for considering this work.

Sincerely,

Zhu Zilong
ORCID: 0000-0002-6955-0903
Department of Neurology, Tianjin Huanhu Hospital, Tianjin, China
zhuzilong1976@gmail.com

---

## Optional block: suggested reviewers

Insert into the submission form if the journal requests reviewers. Rationale is given for each candidate; institutional affiliations and addresses must be verified before submission (not asserted here).

| Candidate | Basis for expertise (verified from the cited literature) | Rationale |
|---|---|---|
| Dr. Shixiang Wu / Dr. Wei Ge | Authors of the systematic evaluation of eight in-silico perturbation methods across four datasets (2026 preprint, doi:10.64898/2026.08.11.744106) | Independent of the tool examined here; positioned to judge the benchmarking context `[affiliation/e-mail to be completed]` |
| Dr. Ru Qiu / Dr. Max Mingqian Zhao | Authors of a confound-diagnostic toolkit for in-silico perturbation (2026 preprint, doi:10.64898/2026.08.04.732812) | Independent; works specifically on confounders of perturbation outputs `[affiliation/e-mail to be completed]` |
| Dr. Kenji Kamimoto | First author of CellOracle (*Nature* 2023; doi:10.1038/s41586-022-05688-9) | Expert on the alternative, multi-hop design used here as the contrast; note this reviewer is an author of a compared tool `[affiliation/e-mail to be completed]` |

Reviewers who authored `scTenifoldKnk` are deliberately **not** suggested, to keep the assessment independent; the editor may of course choose otherwise.

---

## Internal notes (not for submission)

* Letter body: 448 words (within the 300–450 range for an Original Article submission). If the journal imposes a shorter limit, cut paragraph 3 to its first and third sentences (~330 words).
* **Article type: Original Article** (author decision, 2026-09-14). When submitting, match the exact label used in the journal's submission system; if the closest label is "Research Article" or "Original Research", change the phrase "as an Original Article" in paragraph 1 only.
* The letter does not state or imply acceptance likelihood, acceptance rates or editorial preferences.
* If the author list grows, replace "the author" with "all authors" in the Declarations paragraph and add a contributions sentence.
