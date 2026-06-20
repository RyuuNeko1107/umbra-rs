# VSOP87 原データ provenance（ISSUE-033 / `docs/data-sources.md` §2.1, §0/§5）

太陽地心位置のための地球日心黄道座標。**版D**＝黄道・平均分点 of date・球面（L,B,R）。
**版A**＝黄道・分点 **J2000**・直交（X,Y,Z）。版A は M10 太陽フレーム修正（of-date 歳差不整合の解消・
ISSUE-013/035）で採用＝月（ELP, J2000 黄道）と同じ J2000 固定フレーム経路に揃える。

## 取得元（一次配布）

IMCCE（旧 Bureau des Longitudes）/ Bretagnon & Francou (1988), *A&A* 202, 309。

| ファイル | 内容 | URL |
|---|---|---|
| `VSOP87D.ear` | 地球の VSOP87 版D 系列（L, B, R・黄道 of date 球面） | ftp://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87D.ear |
| `VSOP87A.ear` | 地球の VSOP87 版A 系列（X, Y, Z・黄道 **J2000** 直交） | ftp://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87A.ear |
| `vsop87.chk` | 全版・全天体の検証値（評価器 ISSUE-013 の照合用） | ftp://ftp.imcce.fr/pub/ephem/planets/vsop87/vsop87.chk |

- 取得日: `VSOP87D.ear`/`vsop87.chk` = 2026-06-15、`VSOP87A.ear` = 2026-06-20（UTC）。一次配布元 IMCCE。
  各ファイルは公開ミラー（github `ctdk/vsop87` / CDS VizieR `VI/81`）とバイト同一であることを確認
  （`VSOP87D.ear` 324786 B / `VSOP87A.ear` 472948 B、完全性クロスチェック）。
- 取得は開発時の明示手順（本記録）による一回限り。**ビルド時/実行時のネットワーク取得は行わない**
  （`docs/accuracy.md` §5）。`xtask generate-coefficients --dataset vsop87` は本ローカルファイルのみを入力とする。

## 元データ checksum（SHA-256）

```
VSOP87D.ear  8b160c859136d467f2be7fc29efa8a9652e95516dfbde00e4c739d7ddc90ca91
VSOP87A.ear  69d0b4c7525f094a03099e64558321eb64f2402a478472ee537239bcc59b7cb6
vsop87.chk   f8fa52449262be05a22a96840c1acbad0b35c8999e00b5c0477ba8a91a67a51a
```

## 版A・地球系列の明示（XYZ・J2000）

`VSOP87A.ear` ヘッダ各行に `VSOP87 VERSION A1   EARTH   VARIABLE n (XYZ)  *T**α  ... HELIOCENTRIC
DYNAMICAL ECLIPTIC AND EQUINOX J2000`。変数 n∈{1=X,2=Y,3=Z}、α∈{0..5}、18 セクション。全変数 AU。
評価式・項行形式（末尾3値 = A,B,C）は版D と同一だが、結果は J2000 黄道**直交**で球面変換不要。
J2000 黄道（VSOP87 力学的分点）→ ICRS は固定回転（VSOP87→FK5/ICRS）で行う（of-date 歳差を通さない）。

## ライセンス・帰属

VSOP87 係数は IMCCE 由来の科学データ（明示的 SPDX ライセンスなし）。CDS/IMCCE がオープン配布し、
係数は事実データとして帰属付きで多数の許諾型 OSS が数十年同梱（de facto）。本 crate は原一次データから
自前生成する（`docs/data-sources.md` §0/§6）。Bretagnon & Francou (1988) を引用。**GPL 再実装は取り込まない**。

## 版D・地球系列の明示（B4(c) 必須チェック）

ファイルヘッダに `VSOP87 VERSION D4    EARTH` と明記。**版D（黄道 of date 球面 heliocentric）かつ
body=EARTH（EMB ＝地球–月重心ではない）**であることを生成時に検査する。EMB 取り違えは月軌道による
**6.4″/月オーダーの系統誤差**になり、太陽地心＝地球日心の符号反転テストでは検出できない（ISSUE-033/013）。

## 系列の構造（パース仕様の正本）

- ファイルは「セクション見出し行」と「項行」の連続。
  - 見出し: `... VARIABLE v (LBR)  *T**α  N TERMS ...`。v∈{1=L 黄経, 2=B 黄緯, 3=R 動径}、α∈{0..5}、N=項数。
  - 項行: `<id> <rank> <12 整数乗数>  S  K  A  B  C`。**末尾 3 値 = A, B, C**。S,K は別表現（評価不要）。
- 評価式（標準 VSOP87）: 各変数 `s = Σ_{α} T^α · Σ_k A_{α,k}·cos(B_{α,k} + C_{α,k}·T)`。
  T = ユリウス千年 from J2000 TDB = `(JD_TDB − 2451545.0)/365250`。A は L,B で無次元(rad)・R で AU、
  B は rad、C は rad/千年。評価は ISSUE-013（本 Issue は取り込み・packed 化まで）。
- 地球の項数（VSOP87D.ear ヘッダ宣言値）:
  - L: T0=**559**, T1=**341**, T2=**142**, T3=**22**, T4=**11**, T5=**5**
  - B: T0=**184**, T1=**99**, T2=**49**, T3=**11**, T4=**5**（T5 なし）
  - R: T0=**526**, T1=**292**, T2=**139**, T3=**27**, T4=**10**, T5=**3**
  - 合計 **2425 項**。
