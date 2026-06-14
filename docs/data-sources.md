# データ出典とライセンス (data-sources)

`umbra-rs` が使う外部データ（天体暦係数・地球姿勢・閏秒・検証オラクル）の**出典・形式・対応年代・ライセンス・取り込み方針**を定める。
方針は accuracy.md / conventions.md と整合。**ライブラリ本体は実行時にネットワークしない**（取得・生成は xtask に隔離）。

> 状態: ドラフト（Milestone 0）。**OSS公開前にライセンス最終確認（plan §26）が必須**。本書 §6 の未解決事項を潰すまで「ライセンス確認済み」としない。
> 注意: 本書のライセンス記述は調査時点（2026-06）の理解。法的助言ではない。crate へ同梱する前に各一次配布元の最新規約を確認すること。

---

## 0. 最重要: ライセンス上の主要リスク

1. **解析天体暦の係数（VSOP87 / ELP/MPP02）は明示的なOSSライセンスを持たない科学データ**である（IMCCE / 旧 Bureau des Longitudes 由来）。広く再配布・再実装されているが、許諾条件が文書化されていない。
2. **既存の再実装には GPL のものがある**（例: Yuk Tung Liu の ELP/MPP02 C++ 実装と生成係数は **GPL-3.0**）。許諾型（MIT/Apache-2.0 想定）の本crateに**GPLコードや GPL派生係数を取り込んではならない**。
3. → 方針: **一次配布元（IMCCE）の原データを取り込み、係数自体は事実として扱い**、自前の生成パイプライン（xtask）で packed 形式へ変換し、**出典明記＋checksum＋NOTICE** を付す。GPL実装は参照（理解・検証）に留め、コード/データを移植しない。生成データは独立クレート/サブモジュールに隔離し、ライセンス区分を明確にする。

---

## 1. サマリ

| データ | 用途 | プロファイル | 同梱 | 対応年代 | ライセンス状況 |
|---|---|---|---|---|---|
| VSOP87(D) | 太陽（地球公転）位置 | Standard | する（生成係数） | 広域 | 科学データ・明示OSSなし → §2.1 |
| ELP/MPP02 | 月位置 | Standard | する（生成係数） | 広域 | 同上・GPL再実装に注意 → §2.2 |
| JPL DE440/DE441 | Reference暦・差分オラクル | Reference(`jpl`) | しない（任意DL） | DE440:1550–2650 / DE441:-13200–+17191 | 米政府由来・要規約確認 → §2.3 |
| IERS EOP C04 | UT1−UTC・極運動・ΔT | Standard | する（versioned data） | 1962–現在＋予測 | IERS公開データ → §3 |
| 閏秒 (TAI−UTC) | UTC↔TAI | 全 | する | 1972– | IERS/IAU公開 → §3 |
| Espenak–Meeus ΔT | 長期ΔT外挿 | 全 | する（多項式） | 歴史～将来 | NASA TP・公開式 → §3 |
| NASA 5千年日食カタログ | 検証オラクル(第二義) | テスト | fixtures のみ | -1999–+3000 | NASA(公的)・出典明記 → §4 |
| USNO / 各国予報値 | 検証オラクル(第二義) | テスト | fixtures のみ | 公開分 | 出典明記 → §4 |

「同梱」= crate に生成済みデータを含めるか。JPL DEと巨大オラクルは含めない（plan §8）。

---

## 2. 天体暦

### 2.1 VSOP87（太陽 / 惑星）
- 出典: Bretagnon & Francou (1988), *A&A* 202, 309。配布: IMCCE（旧 Bureau des Longitudes）FTP / CDS VizieR。
- 採用版: **VSOP87D**（黄道・平均分点 of date、球面）。太陽地心位置は地球の日心位置を反転して得る。
- 形式: 各変数（L, B, R）の周期項テーブル（振幅・位相・振動数）。
- 対応年代: 数千年規模で有効（現代付近で最良）。本プロジェクト v1.0 範囲 1900–2100 は十分内側。
- 精度方針: DE 差分で**残差 0.05″ 級**まで項を採用（accuracy.md §2.4 / §3.3）。
- ライセンス: 明示OSSなし。**科学データとして取り込み、Bretagnon & Francou を引用**。原ファイルから自前生成し provenance + checksum を記録。

### 2.2 ELP/MPP02（月）
- 出典: Chapront & Francou (2002)「The lunar theory ELP revisited」, IMCCE。MPP02 は LLR 版と DE405/406 fit 版の2パラメータ系列を提供。
- 採用版: **DE fit 版**（本プロジェクトの Reference=JPL DE と整合させるため）。
- 形式: 主問題（main problem）級数＋惑星/地球摂動級数。経度・緯度・距離。
- 精度方針: DE 差分で**残差 0.1″ 級**まで項を採用（accuracy.md §2.4）。打切り次数・達成残差を `EphemerisMetadata` に記録。
- ライセンス: **重要** — 一次データは明示OSSなし（科学データ扱い）。**Yuk Tung Liu の再実装/係数は GPL-3.0 のため本crateに取り込まない**（参照・検証のみ可）。原 IMCCE データから自前生成し、Chapront & Francou を引用。
- 速度: 解析微分 or 対称差分でバックエンドが供給（architecture §4）。

### 2.3 JPL DE（Reference・差分オラクル）
- 出典: Park, Folkner, Williams, Boggs (2021), *AJ* 161:105, DOI:10.3847/1538-3881/abd414。
- 配布: NAIF SPK（`.bsp`）`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/`。ASCII 版は JPL SSD。
- 対応年代: **DE440 = 1550–2650**（v1.0 範囲を包含、既定 Reference）。DE441 = -13200–+17191（将来の長期拡張用）。
- 取り込み: feature `jpl` でのみ有効。**crate に同梱しない**。利用者が任意取得し、本crateは SPK reader（または既存純Rust reader の評価）で読む。reader実装の選定は Milestone 10。
- ライセンス: 米政府/Caltech-JPL。データ利用は広く行われているが**OSS同梱の可否は規約確認が必要**（§6）。同梱しない方針によりリスクを回避。引用は上記論文。

---

## 3. 時刻・地球姿勢

### 3.1 IERS EOP C04（UT1−UTC・極運動）
- 出典/配布: IERS Earth Orientation Center, `https://datacenter.iers.org/`（VLBI/GNSS/SLR/DORIS 統合）。
- 系列: **現行運用 C04 系列を使用**（旧 `EOP 14 C04`、ITRF2020対応の新 `EOP 20 C04` へ移行中。採用版を metadata に固定）。1日間隔、1962–現在、+短期予測。
- 用途: UT1−UTC（ERA・恒星時）、極運動（ITRS 変換・中心線位置）。accuracy.md §2.3。
- ライセンス: IERS 公開科学データ。**出典・取得日・系列バージョンを記録**して同梱（versioned + checksum）。

### 3.2 閏秒（TAI−UTC）
- 出典: IERS Bulletin C / IANA tz の `leap-seconds.list`。1972– の積算閏秒。
- 用途: UTC↔TAI↔TT。`data/leap-seconds/` に versioned 同梱。期限切れ検知（valid_to）を持つ。

### 3.3 ΔT（= TT − UT1）
- 履歴/近傍: EOP C04 から導出（高精度）。
- 長期/外挿: **Espenak–Meeus 多項式**（NASA、1620以前および将来外挿）。短期将来は IERS Bulletin A 予測。
- **不確実性帯**を `CalculationMetadata.delta_t_uncertainty_seconds` に出力（将来 UTC 律速、accuracy.md §0/§2.3）。
- ライセンス: Espenak–Meeus 式は NASA 公開（出典明記）。

---

## 4. 検証オラクル（fixtures のみ・本体非依存）

accuracy.md §3.1 の階層に従い、**第一義は JPL DE 差分**（§2.3）。以下は第二義（整合チェック）。

### 4.1 NASA 5千年日食カタログ（Espenak & Meeus）
- 出典: NASA「Five Millennium Catalog of Solar Eclipses」`eclipse.gsfc.nasa.gov`、NASA/TP-2006-214141 等。-1999–+3000。
- 用途: 種別・最大食時刻・gamma・食分・ベッセル要素の整合チェック。
- 注意: NASA は固有の **ΔT・k 慣習**（Espenak 2値 k、conventions §9）を用いる。**慣習を揃えた上で比較**し、系統差を accuracy.md に記録。絶対基準にしない。
- ライセンス: NASA（米政府）成果。**出典明記**で fixtures に取り込み。サイト掲載物の体裁は引用元を尊重。

### 4.2 USNO / 各国機関の予報値
- 用途: 局地予報・接触時刻のスポット整合。
- 取り込み: 値を fixtures に転記し**出典・取得日を併記**。商用配布物・著作性のある表現は転記しない（数値事実のみ）。

---

## 5. データ取り込みパイプラインと完全性

```
data/coefficient-source/   ← 一次配布元の原データ（VSOP87 / ELP-MPP02 等。再配布可否は出典に従う）
   ↓ cargo xtask generate-coefficients   （打切り・packed化・provenance付与）
generated/                 ← crate 組込み用バイナリ/Rust
   ↓ cargo xtask verify-generated         （CIで生成済みとの差分検査）
crate へ組込み
```

- 各データセットに `DataSetMetadata { name, version, source, license, valid_from, valid_to, checksum }`（architecture §11）。
- 再現性: 原データ→生成物が決定的。CI で checksum 差分を検出。
- 隔離: 生成係数は**ライセンス区分が明確なクレート/モジュール**に置き、NOTICE・引用を同梱。GPL派生物を混入させない（§0）。

---

## 6. ライセンス調査結果（2026-06、一次情報による）

> まとめ: **コードは MIT/Apache 維持で問題なし。同梱データは「オープン配布 + 帰属表示」の科学データとして再配布可（法的リスクは低い）。ただし VSOP87/ELP には crisp な SPDX ライセンスが付かない点と、GPL 再実装の非混入が要注意。** 詳細は各項。

### 6.1 各データ源の判定

- **JPL DE440/441 — 明確にパブリックドメイン（問題なし）。** NAIF 公式: 「All of the data and tools generated by NAIF are free and in the public domain.」 再配布・改変自由。本プロジェクトは Reference 用で同梱しない方針だが、同梱しても可。出典: Park et al. (2021) を引用。
- **VSOP87（CDS VizieR VI/81 / Bretagnon & Francou 1988）— 再配布可（帰属必須）。低リスク。** CDS の法務（cds.unistra.fr/legals）はカタログを **Open Licence / ODbL / CC-BY 等のオープンライセンス**で配布し DOI/出典明記を求める。ただし VI/81 個別ページに明示ライセンス表示はなく、データファイル自体にもライセンス文言はない。係数（フーリエ振幅・振動数）は**科学的事実データ**で、米法では事実は非著作物・EU は sui generis DB 権があるが CDS がオープン配布。**多数の許諾型（MIT/Apache/BSD）OSS が数十年 VSOP87 を同梱**（de facto）。→ Bretagnon & Francou (1988) + CDS を帰属し、データとして同梱可。
- **ELP/MPP02（IMCCE FTP / Chapront & Francou 2002, 2003）— 同上（帰属必須・低リスク）。** 原データ（ELP_MAIN/PERT.S*）は IMCCE FTP 配布、ファイルに明示ライセンスなし＝科学データ扱い。**重要: ytliu0/ElpMpp02（GPL-3.0）と MarcvdSluys/ELP-MPP02 等の再実装/派生係数は取り込まない**（参照のみ）。CALCEPH ソフトは CeCILL だが本件はデータの話で無関係。→ **原 IMCCE データから自前生成**し Chapront & Francou を引用。
- **IERS EOP / 閏秒 — オープンアクセス（帰属）。** IERS はオープンデータ方針（re3data: open）。明示ライセンス文は薄いが、EOP は普遍的に再配布されており、閏秒（TAI−UTC）は IANA tz・各 OS に普遍同梱の事実データ。→ 引用付きで同梱可。
- **Espenak–Meeus ΔT / NASA 5千年カタログ — NASA(米政府)成果・公開式。** ΔT 多項式は公開数式（既に実装、ライセンス問題なし）。カタログは出典明記で fixtures 転記可。

### 6.2 結論と運用方針

1. **コードライセンス**: MIT OR Apache-2.0（暫定確定可）。
2. **同梱データ**: 「オープン配布の科学データ + 帰属」。crisp な SPDX ライセンスが付かないため、`LICENSE-*` ではなく **NOTICE / data ディレクトリの出典表記**で各データセットの出典・引用・配布元を明示する（`DataSetMetadata.license` は "scientific data, redistributed with attribution; see NOTICE" 等）。
3. **GPL 非混入（必須）**: VSOP/ELP の GPL 再実装・派生係数を取り込まない。**原一次配布元から自前パイプライン生成**。`cargo-deny` のライセンス allow-list で機械チェック（コード依存。データは別管理）。
4. **JPL DE は非同梱**（feature `jpl` で任意DL、Reference 専用）。
5. **EU sui generis DB 権の残リスク（低）**: 念のため最大安全策として「公開理論からの自前再生成（clean-room）＋出典明記」で運用すれば、単なるファイルコピーより安全。実務上は公開係数の同梱が業界標準。

### 6.3 残チェック（実装時）

- [x] VSOP87 / ELP/MPP02 再配布可否 → **可（帰属付き・低リスク。GPL再実装は不可）**
- [x] JPL DE → **PD・非同梱方針で可**
- [x] IERS EOP / 閏秒 → **オープン・帰属付きで可**（valid_to 超過時 `Missing*Data` は実装済み）
- [ ] NOTICE ファイル作成（全データ源の出典・引用・配布元・各データの扱いを明記）
- [ ] `cargo-deny` allow-list 運用（コード依存）＋ データは NOTICE 管理（ISSUE-046, CI 済み）
- [ ] 係数生成パイプライン（ISSUE-033/034）で原一次データのみ使用・provenance/checksum 記録

---

## 参照
- VSOP87: Bretagnon & Francou (1988), A&A 202, 309。IMCCE / CDS VizieR。
- ELP/MPP02: Chapront & Francou (2002), IMCCE。参考実装 ytliu0/ElpMpp02（**GPL-3.0、取り込み不可**）。
- JPL DE440/441: Park et al. (2021), AJ 161:105。NAIF。
- IERS EOP C04 / Bulletin A/C: `iers.org` / `datacenter.iers.org`。
- NASA 日食カタログ: Espenak & Meeus, `eclipse.gsfc.nasa.gov`。
