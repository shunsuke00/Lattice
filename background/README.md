これは単一場インフレーションの背景場の発展を調べるためのプログラムで
"Simulating Inflationary Universe"を再現するように作った。

[理論的セットアップ]
・換算プランク質量M_pl=1/sqrt(8*pi*G)=1とする単位系で考える
・スケール因子の初期値は1に取る

[インフラトンポテンシャル]
・調和振動子(カオスインフレーション)
・stepポテンシャル
の2種類を用意。model.hppで実装されている。

[リスケーリング]
・phi_pr = phi [M_pl]
・x_pr   = m x
・dt_pr  = m a^s dt

[時間積分]
・4次のルンゲクッタ法
・オイラー法
の2種類を用意。evolution.hppで実装されている。
どの積分法を用いるかはparameter.hppのMETHODマクロで選択する。

[パラメータの選択]
・dtは10^-4オーダーで取る
・スケール因子が1000まで見る
    カオスインフレーション + 宇宙時間 + s=1 --- dt=1.0e-4,total_step=1.8e+6
    カオスインフレーション + 共形時間 + s=1 --- dt=1.0e-4,total_step=1.2e+4
    step potential + 宇宙時間 --- 
    step potential + 共形時間 + s=1 --- dt=1.0e-4,total_step=

[出力]
・ハッブルパラメータは宇宙時間でのものを出力
・場の微分は宇宙時間での微分を出力

[グラフの描画]
graph.pyのpythonファイルをrunさせると自動でグラフを出力してくれる
必要なものはmatplotlib

[インフラトンポテンシャル]
potential.pyはインフラトンポテンシャルの形を確認するための簡易的なpythonファイル
独立で実行可能

[初期条件]
論文で与えられている初期条件を求めるためのプログラムを別個で"background_initial"ディレクトリ
に作ったが実際は必要なく、このプログラムで
f=15.0
df=0.0
として走らせてやれば良い。

[未遂]
・プログラム上の時間間隔を a^s-1 dt_0 に取る方法はまだ実装していない
->実装していないうちはs=1で固定しておくべき
・リープフロッグ法の実装はまだ
・axion-U(1)の背景場プログラムはAがdeltaphiから生成されるので不可能