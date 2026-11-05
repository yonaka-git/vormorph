# 参考文献

この文書は、vormorph の設計と実装で参照する文献の一覧です。

- 書誌情報は 2026-11-05 に、Crossref、J-STAGE、doi.org、出版社・公式サイトで確認しました。
- オープンアクセスの文献は `references/` に保存しています。著作権のためリポジトリには含めず、`.gitignore` で除外しています。
- 入手状況の欄の意味は次のとおりです。

| 表記 | 意味 |
|---|---|
| 取得済み | 正規の配布元から取得し、`references/` に保存済み |
| 要手動入手（有料） | 出版社の有料コンテンツ |
| 要手動入手（無料公開） | 出版社が無料公開しているが、自動では取得できなかった |
| 書誌のみ | オンラインの全文がない（会議論文・学位論文・書籍など） |
| URL のみ | ソフトウェアのリポジトリなど |

## 半陰解法・wet/dry

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Casulli, V., Walters, R. A. | 2000 | An unstructured grid, three-dimensional model based on the shallow water equations | Int. J. Numer. Meth. Fluids, 32(3), 331–348 | [10.1002/(SICI)1097-0363(20000215)32:3<331::AID-FLD941>3.0.CO;2-C](https://doi.org/10.1002/(SICI)1097-0363(20000215)32:3%3C331::AID-FLD941%3E3.0.CO;2-C) | 要手動入手（有料） |
| Casulli, V. | 2009 | A high-resolution wetting and drying algorithm for free-surface hydrodynamics | Int. J. Numer. Meth. Fluids, 60(4), 391–408 | [10.1002/fld.1896](https://doi.org/10.1002/fld.1896) | 要手動入手（有料） |
| Casulli, V., Zanolli, P. | 2012 | Iterative solutions of mildly nonlinear systems | J. Comput. Appl. Math., 236(16), 3937–3947 | [10.1016/j.cam.2012.02.042](https://doi.org/10.1016/j.cam.2012.02.042) | 要手動入手（無料公開。Elsevier open archive。自動取得は権限チェックで止まったため未取得） |

## 急変流・保存型スキーム

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Stelling, G. S., Duinmeijer, S. P. A. | 2003 | A staggered conservative scheme for every Froude number in rapidly varied shallow water flows | Int. J. Numer. Meth. Fluids, 43(12), 1329–1354 | [10.1002/fld.537](https://doi.org/10.1002/fld.537) | 要手動入手（有料） |
| Kramer, S. C., Stelling, G. S. | 2008 | A conservative unstructured scheme for rapidly varied flows | Int. J. Numer. Meth. Fluids, 58(2), 183–212 | [10.1002/fld.1722](https://doi.org/10.1002/fld.1722) | 要手動入手（有料） |
| Kernkamp, H. W. J., van Dam, A., Stelling, G. S., de Goede, E. D. | 2011 | Efficient scheme for the shallow water equations on unstructured grids with application to the Continental Shelf | Ocean Dynamics, 61(8), 1175–1188 | [10.1007/s10236-011-0423-6](https://doi.org/10.1007/s10236-011-0423-6) | 要手動入手（有料。TU Delft Repository に投稿版の記録あり） |

## 非構造 C 格子・再構成

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Perot, B. | 2000 | Conservation properties of unstructured staggered mesh schemes | J. Comput. Phys., 159(1), 58–89 | [10.1006/jcph.2000.6424](https://doi.org/10.1006/jcph.2000.6424) | 要手動入手（有料） |
| Ringler, T. D., Thuburn, J., Klemp, J. B., Skamarock, W. C. | 2010 | A unified approach to energy conservation and potential vorticity dynamics for arbitrarily-structured C-grids | J. Comput. Phys., 229(9), 3065–3090 | [10.1016/j.jcp.2009.12.007](https://doi.org/10.1016/j.jcp.2009.12.007) | 要手動入手（有料） |

## 格子生成

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Engwirda, D. | 2017 | JIGSAW-GEO (1.0): locally orthogonal staggered unstructured grid generation for general circulation modelling on the sphere | Geosci. Model Dev., 10(6), 2117–2140 | [10.5194/gmd-10-2117-2017](https://doi.org/10.5194/gmd-10-2117-2017) | 取得済み（`engwirda2017_jigsaw-geo.pdf`、CC BY 3.0） |
| Engwirda, D. | — | JIGSAW（C++ の格子生成ライブラリ） | GitHub | <https://github.com/dengwirda/jigsaw> | URL のみ（独自ライセンス。商用システムへの組み込みには作者との取り決めが必要） |
| Engwirda, D. | — | jigsaw-python（JIGSAW の Python インターフェース） | GitHub | <https://github.com/dengwirda/jigsaw-python> | URL のみ（同上） |

## 流砂

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Meyer-Peter, E., Müller, R. | 1948 | Formulas for bed-load transport | Proc. 2nd Meeting of the International Association for Hydraulic Structures Research (IAHSR), Stockholm, pp. 39–64 | DOI なし | 書誌のみ |
| 芦田和男, 道上正規 | 1972 | 移動床流れの抵抗と掃流砂量に関する基礎的研究 | 土木学会論文報告集, 第206号, 59–69 | [10.2208/jscej1969.1972.206_59](https://doi.org/10.2208/jscej1969.1972.206_59) | 取得済み（`ashida1972_bedload.pdf`） |
| 長谷川和義 | 1984 | 沖積蛇行の平面および河床形状と流れに関する水理学的研究 | 北海道大学 博士論文 | [NDL サーチ](https://ndlsearch.ndl.go.jp/books/R100000002-I000000198934) | 書誌のみ（紙媒体。勾配補正の式の原著として Nays2DH マニュアルなどが引用） |
| 長谷川和義 | 1983 | 蛇曲流路における流れと平衡底面形状に関する研究 | 土木学会論文報告集, 第338号, 105–114 | [10.2208/jscej1969.1983.338_105](https://doi.org/10.2208/jscej1969.1983.338_105) | 取得済み（`hasegawa1983_meanderbed.pdf`。勾配補正の係数の式を学会誌で示した論文） |
| 長谷川和義 | 1981 | 非平衡性を考慮した側岸浸食量式に関する研究 | 土木学会論文報告集, 第316号, 37–50 | [10.2208/jscej1969.1981.316_37](https://doi.org/10.2208/jscej1969.1981.316_37) | 取得済み（`hasegawa1981_bankerosion.pdf`。斜面上の流砂の導出） |
| Garcia, M., Parker, G. | 1991 | Entrainment of bed sediment into suspension | J. Hydraul. Eng., 117(4), 414–435 | [10.1061/(ASCE)0733-9429(1991)117:4(414)](https://doi.org/10.1061/(ASCE)0733-9429(1991)117:4(414)) | 要手動入手（有料） |
| Itakura, T., Kishi, T.（板倉忠興, 岸力） | 1980 | Open channel flow with suspended sediments | J. Hydraul. Div., ASCE, 106(HY8), 1325–1343 | [10.1061/JYCEAJ.0005483](https://doi.org/10.1061/JYCEAJ.0005483) | 要手動入手（有料。板倉・岸の浮上量式の原著） |
| 板倉忠興 | 1984 | 河川における乱流拡散現象に関する研究 | 北海道開発局 土木試験所報告, 第83号 | [NDL デジタルコレクション](https://dl.ndl.go.jp/pid/8706675) | 書誌のみ（NDL のビューアで閲覧可。浮上量式を和文で扱った文献） |
| Rubey, W. W. | 1933 | Settling velocity of gravel, sand, and silt particles | American Journal of Science, Ser. 5, 25(148), 325–338 | [10.2475/ajs.s5-25.148.325](https://doi.org/10.2475/ajs.s5-25.148.325) | 取得済み（`rubey1933_settling-velocity.pdf`） |
| 平野宗夫 | 1971 | Armoring をともなう河床低下について | 土木学会論文報告集, 第195号, 55–65 | [10.2208/jscej1969.1971.195_55](https://doi.org/10.2208/jscej1969.1971.195_55) | 取得済み（`hirano1971_armoring.pdf`。交換層モデルの原著） |
| 平野宗夫 | 1972 | 混合砂れき河床の変動と平衡河床に関する研究 | 土木学会論文報告集, 第207号, 51–60 | [10.2208/jscej1969.1972.207_51](https://doi.org/10.2208/jscej1969.1972.207_51) | 書誌のみ記録（J-STAGE で公開。1971 年の続報） |

## 実装・マニュアル（公式配布）

| 著者 | 年 | 題名 | 発行 | URL | 入手状況 |
|---|---|---|---|---|---|
| iRIC Project | 2023 | Nays2DH ソルバーマニュアル（日本語版） | iRIC Project | [配布ページ](https://i-ric.org/solvers/nays2dh/) | 取得済み（`iric2023_nays2dh_manual_ja.pdf`。ファイル名上 v4、2023-06-13 版） |
| Deltares | 2026 | D-Flow Flexible Mesh Technical Reference Manual | Deltares | [配布ページ](https://oss.deltares.nl/web/delft3dfm/downloads) | 取得済み（`deltares2026_dflowfm-tech-ref-manual.pdf`。Version 2027.01, Revision 81397, 2026-11-02 の DRAFT 版。随時更新される） |
| MPAS Project | 2015 | MPAS Mesh Specification, Version 1.0 | MPAS Project | <https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf> | 取得済み（`mpas2015_mesh-spec.pdf`） |

## 検証ケースで使う文献

| 著者 | 年 | 題名 | 誌名・巻(号)・頁 | DOI / URL | 入手状況 |
|---|---|---|---|---|---|
| Thacker, W. C. | 1981 | Some exact solutions to the nonlinear shallow-water wave equations | J. Fluid Mech., 107, 499–508 | [10.1017/S0022112081001882](https://doi.org/10.1017/S0022112081001882) | 要手動入手（有料） |
| Ritter, A. | 1892 | Die Fortpflanzung der Wasserwellen | Zeitschrift des Vereines Deutscher Ingenieure, 36 | DOI なし | 書誌のみ（頁は文献により 947–954 または 934–954 と記載が分かれ、原本は未確認） |
| Stoker, J. J. | 1957 | Water Waves: The Mathematical Theory with Applications | Interscience（Wiley Classics Library 版 1992） | [10.1002/9781118033159](https://doi.org/10.1002/9781118033159) | 要手動入手（書籍） |

## 書誌情報の確認メモ

当初の文献リスト（記憶にもとづく書誌情報）から、次の点を修正・特定しました。

- **Meyer-Peter & Müller (1948)**：会議の名称は「2nd IAHR Congress」ではなく、IAHR の前身である IAHSR の第2回会議（Stockholm, 1948年6月7〜9日）です。頁は pp. 39–64 です（Wong & Parker, 2006 の参考文献欄で確認）。
- **芦田・道上 (1972)**：共著者の表記は J-STAGE では「道上 正規」です。誌名は「土木学会論文報告集」です（Nays2DH マニュアルの参考文献欄は「土木学会論文集」としています）。
- **長谷川の式**：Nays2DH マニュアルと Shimizu & Itakura (1989) が引用する原著は、長谷川 (1984) の北海道大学博士論文です。学会誌で勾配補正の係数 $\gamma = \sqrt{\tau_{\ast c}/(\mu_s\mu_k\tau_\ast)}$ が示されているのは長谷川 (1983) で、その導出は長谷川 (1981) を参照しています。Nays2DH では、この補正を河床近傍の流速の向きと組み合わせて使っています。
- **板倉・岸の浮上量式**：原著は英文の Itakura & Kishi (1980) です。和文の土木学会論文としての原著は見つかりませんでした。和文で式を扱った文献としては、板倉 (1984) があります。
- **平野 (1971)**：交換層モデルの原著は「Armoring をともなう河床低下について」（第195号）です。「拡幅を伴う流路変動について」は 1973 年（第210号）の別の論文で、交換層の原著ではありません。
- **Rubey (1933)**：沈降速度式の原著を特定しました（American Journal of Science）。
- **D-Flow FM のマニュアル**：公式の配布ファイルは随時更新される DRAFT 版です。引用するときは版とリビジョンを併記します。

## 追加候補

検証ケース V08〜V11 で比較対象にする実験・解析は、実装時に選定してこの文書に追加します。
