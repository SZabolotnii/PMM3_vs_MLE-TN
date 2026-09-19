# TODO по репозиторію після скорочення статті до 10 сторінок

Контекст: `submission/SZabolotnii_ITEST_paper_10p.docx` зроблено з
`submission/SZabolotnii_ITEST_paper+++ref.docx` вручну (правки XML), а не з R-пайплайну.
Перелік самих правок — у `submission/CHANGES_14p_to_10p.md`.

Стан на 2026-09-19: зроблено все, що робиться в репозиторії.
Лишилось те, що потребує рук автора: Zenodo, Word і **правки тексту статті (розд. 0 нижче)**.

---

## 0. Помилки в тексті статті, знайдені при перевірці відтворюваності (блокує подачу)

Чистий прогін `scripts/run_all.R` (без кешу, seed 20260312) показав, що
`T6_misspecification.csv` і `T7_bimodal_test.csv` у первинному коміті
(eb287cf, 2026-03-13) були пораховані заглушкою PMM3, яка повертала OLS:
ARE(PMM3) = 1 рівно в усіх рядках T6, `mean_iter = 0` в усіх рядках T7.
Колонки MLE-TN у T6 відтворюються точно, тобто дані ті самі — різниться лише PMM3.
T1 і T2 відтворюються біт-у-біт (T2 з точністю 1e-17).
Обидва CSV замінено на вихід пайплайну; `paper/tables/*.tex` перегенеровано.

Нова Табл. 2 (робастність, n = 100, M = 500):

| True DGP | γ₄ | ARE(MLE-TN) | ARE(PMM3) було | ARE(PMM3) стало |
|---|---|---|---|---|
| Uniform | −1.2 | 1.655 | 1.000 | **2.685** |
| Triangular | −0.6 | 0.969 | 1.000 | **1.061** |
| Logistic | +1.2 | 0.994 | 1.000 | **1.023** |
| Student's t10 | +1.0 | 0.982 | 1.000 | **1.047** |

Нові значення узгоджуються з теорією g₃ = 1 − γ₄²/(6 + 9γ₄ + γ₆): для t10 теоретичне
ARE ≈ 1.04, для logistic ≈ 1.07, для uniform ≈ 3.3 (асимптотично; при n = 100 — 2.7).

Що переписати в `submission/SZabolotnii_ITEST_paper_10p.docx`:
- [ ] **Табл. 2**, колонка ARE(PMM3) — значення вище.
- [ ] **Розд. 3.3, абзац після Табл. 2** — «PMM3 is robust, maintaining ARE ≈ 1.0 for all
      investigated distributions». Насправді PMM3 не гірший за OLS в усіх чотирьох
      випадках, а на Uniform (γ₄ = −1.2) дає ARE 2.68 — краще, ніж MLE-TN (1.655).
      Висновок «PMM3 — безпечний вибір» лишається і навіть посилюється.
- [ ] **Розд. 3.3, останнє речення** про збіжність — «100% of the 500 replications …
      absolute bias never exceeding 0.022». Пайплайн запускає цей тест з M = 200
      (`run_all.R`, Block 6). Нові числа: 100% збіжності, 0% NA, 2.8–3.7 ітерації
      в середньому, |bias| ≤ 0.003.
- [ ] **Внесок (iii) у вступі** — «PMM3 maintains ARE ≈ 1 under model misspecification».
- [ ] **Conclusions, пункт «Robustness»** — те саме «ARE ≈ 1».
- [ ] **Abstract** — «PMM3 is shown to be robust under misspecification» можна лишити,
      але варто додати число (ARE 1.02–2.68).

Інші розбіжності тексту з даними (від пайплайну не залежать):
- [ ] **Межа бімодальності.** Розд. 2.2 («λ > 1/√2 — the distribution becomes bimodal»),
      розд. 2.6 («bimodal for λ > 1/√2, that is for every λ ≥ 1.5») і підпис Рис. 1
      («bimodality boundary λ = 1/√2»). Для TN(λ) = ½N(−λ,1) + ½N(λ,1) межа — **λ = 1**
      (друга похідна густини в нулі: λ² − 1). Саме так рахує `tn_cumulants()`
      і Табл. T1 (λ = 1 — unimodal). Рисунок уже виправлено, текст — ні.
- [ ] **Iris, розд. 3.2**: «g₃ theor = 0.517 which corresponds to a 51.7% reduction».
      Має бути g₃ theor = 0.483 (1 − 0.483 = 51.7%).
- [ ] **Negative example, розд. 3.2**: «four further real datasets with γ₄ ∈ [−0.54, −0.12],
      for which g₃ emp ∈ [0.97, 1.03]». За `real_data_all_candidates.csv` і звітом
      (розд. 9.5) у trees / cars / faithful / iris-setosa / biaxial PMM3 має
      bootstrap g₃ emp = 1.28 / 1.14 / 1.03 / 1.22 / 1.36, тобто ∈ [1.03, 1.36] —
      не краще, а до 36% гірше за OLS; LOO-MSE/OLS ∈ [1.000, 1.013].
      Поріг γ₄ ≲ −0.7 це підтримує ще сильніше, але формулювання «equivalent to OLS»
      у Conclusions (п. 3) варто уточнити: за прогнозом — так, за дисперсією — ні.
- [ ] Iris і real-data CSV не генеруються жодним скриптом у репо (`run_all.R` їх
      не пише). Числа в статті з ними збігаються, але відтворити їх з коду не можна.
      Бажано додати скрипт, що їх рахує (B = 2000, seed 20260312).

---

## 1. Супплемент (блокує подачу)

- [x] `supplement/` зібрано скриптом `scripts/09_build_supplement.R`: 8 CSV
      (T1, T2, T6, T7, iris, три real_data), `README.md` з описом колонок,
      `SESSION_INFO.txt` (версії R і пакетів), `MD5SUMS.txt`
- [ ] залити на Zenodo, отримати DOI, вписати його в `supplement/README.md`
- [ ] вставити DOI у два місця в docx: розд. 3.1 («the full numeric table is provided
      as supplementary material») і в кінці розд. 3.2

---

## 2. Генерація рисунків

- [x] прибрано `title =` і `subtitle =` з усіх `labs()` у `scripts/08_generate_paper_figures.R`
- [x] двопанельний рисунок через patchwork: `paper/figures/Fig1_2_TN_densities_efficiency.png`,
      24×7 см, панелі (a)/(b); окремі Fig1/Fig2 теж генеруються
- [x] межу бімодальності на панелі (b) перенесено з λ = 1/√2 на λ = 1 (див. розд. 0)
- [x] Fig3_g3_convergence, Fig4_ARE_comparison — без заголовків
- [ ] перевставити в docx: Рис. 1 ← `Fig1_2_TN_densities_efficiency.png`,
      Рис. 2 ← `Fig3_g3_convergence.png`, Рис. 3 ← `Fig4_ARE_comparison.png`;
      оновити підпис Рис. 1 (λ = 1)

---

## 3. Таблиці генеруються, а не набираються руками

- [x] `scripts/07_generate_paper_tables.R` тепер пише ще `paper/tables/paper_tables.md`
      і `paper/tables/paper_tables.docx` (через pandoc) з тих самих data frame, що й .tex
- [x] звірено Табл. 1 (iris) статті з CSV — збігається
- [ ] Табл. 2 (робастність) — перенести нові значення з `paper_tables.docx` (розд. 0)

---

## 4. Узгодити набір n

Нічого робити не треба: усі чотири n ∈ {46, 100, 200, 500} видно на рисунках.
Не повертати таблицю Монте-Карло в скороченому вигляді.

---

## 5. Гігієна git

- [x] `.gitignore`: `~$*.docx`, `~$*.xlsx`, `~$*.pptx`, `submission/_preview/`
- [x] з репо прибрано lock-файл Word і `submission/_preview/`
- [x] обидві версії статті + CHANGES закомічені раніше (ac8ed5c)

---

## 6. Перед відправкою

- [x] `Rscript scripts/run_all.R` начисто (без `results/mc_cache`) — ~75 с.
      Заодно виправлено: запуск через `Rscript` шукав `R/config.R` на рівень вище
      і падав; блок рисунків падав, бо `results/figures/` немає у свіжому клоні.
- [ ] після правок розд. 0 відкрити 10p.docx у Word, Ctrl+A → F9 (поля SEQ);
      має вийти Табл. 1–2, Рис. 1–3
- [ ] перевірити обсяг у Word — 10 сторінок пораховано конвертацією LibreOffice,
      на останній сторінці ~35% вільного місця
