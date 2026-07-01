# 🧭 OmicsCompass

オミクスデータから経路を探索する研究支援ツール。RNA-seq などの公開データを取得し、差分解析・パスウェイ解析・ネットワーク解析・上流因子推定を行うことができます。

## 🌐 アクセス URL

https://omics-compass-emrvfcsqsf8wgztwhdgjyd.streamlit.app

ログインが必要です（研究室内で共有する場合は管理者に認証情報を確認してください）。

---

## ✨ 機能一覧

| ページ | 機能 |
|--------|------|
| 📥 Data Fetch | NCBI GEO から RNA-seq・マイクロアレイデータを取得 |
| 🔬 Analysis | 差分発現解析（t検定 + Benjamini-Hochberg 補正）、プローブ→遺伝子シンボル変換、LLM によるグループ自動判定 |
| 📊 Visualization | Volcano Plot、Heatmap、パスウェイ解析、遺伝子ネットワーク解析（STRING DB）、上流転写因子推定 |

---

## 🚀 ローカルでの起動方法

### 必要なもの

- Python 3.12+
- Git
- （オプション）Ollama（ローカル LLM）

### セットアップ

```bash
# リポジトリをクローン
git clone https://github.com/makoisaiah/omics-compass.git
cd omics-compass

# 仮想環境を作成・有効化
python3 -m venv venv
source venv/bin/activate

# パッケージをインストール
pip install -r requirements.txt

# Streamlit Secrets を設定
mkdir -p .streamlit
cp .streamlit/secrets.toml.example .streamlit/secrets.toml
# secrets.toml を編集して認証情報を入力

# アプリを起動
streamlit run app.py
```

### Secrets の設定

`.streamlit/secrets.toml` に以下の形式で記述してください。

```toml
[credentials.usernames.あなたのユーザー名]
name = "表示名"
password = "ハッシュ化されたパスワード"

[cookie]
name = "omics_compass_cookie"
key = "32文字以上のランダムな文字列"
expiry_days = 30
```

パスワードのハッシュ化:

```bash
python3 -c "
import streamlit_authenticator as stauth
credentials = {
    'usernames': {
        'yourusername': {'password': 'yourpassword'}
    }
}
print(stauth.Hasher.hash_passwords(credentials))
"
```

---

## 🔬 使い方

### 1. データ取得

1. `Data Fetch` ページを開く
2. NCBI GEO（https://www.ncbi.nlm.nih.gov/geo/）で目的のデータセットを検索
3. GSE 番号（例: `GSE21393`）を入力して取得

### 2. 差分解析

1. `Analysis` ページを開く
2. 「LLM でグループを自動判定」ボタンでコントロール群・処理群を自動推定
3. 必要に応じて手動で修正
4. 「差分解析を実行」ボタンで解析開始

### 3. 可視化・解析

`Visualization` ページの各タブで解析を実行:

- **Volcano Plot**: 発現変化の概観
- **Heatmap**: 上位遺伝子の発現パターン
- **Pathway Analysis**: KEGG・Reactome・GO・WikiPathways
- **Network**: STRING DB を使った遺伝子相互作用ネットワーク
- **Upstream**: 上流転写因子の推定

---

## 🤖 LLM によるグループ自動判定

ローカルに Ollama が起動していれば自動的にローカル Mistral を使用します。なければ Groq API にフォールバックします。

```bash
# Ollama のインストール
curl -fsSL https://ollama.com/install.sh | sh

# Mistral モデルのダウンロード（約4.1GB）
ollama pull mistral

# 起動確認
ollama run mistral "Hello"
```

---

## 📦 使用パッケージ

| パッケージ | 用途 |
|-----------|------|
| streamlit | Web アプリフレームワーク |
| streamlit-authenticator | ログイン認証 |
| GEOparse | NCBI GEO からデータ取得 |
| pandas / numpy | データ処理 |
| scipy / statsmodels | 統計解析 |
| plotly | インタラクティブ可視化 |
| gseapy | パスウェイ解析・上流因子推定 |
| mygene | プローブ ID → 遺伝子シンボル変換 |
| networkx | ネットワーク解析 |
| requests | STRING DB API アクセス・Ollama 通信 |

---

## 📁 ディレクトリ構成

```
omics-compass/
├── app.py                    # メインアプリ・ログイン
├── requirements.txt          # 依存パッケージ
├── .python-version           # Python バージョン指定（3.12）
├── .gitignore
├── pages/
│   ├── 01_data_fetch.py     # GEO データ取得
│   ├── 02_analysis.py       # 差分解析
│   └── 03_visualization.py  # 可視化・各種解析
├── components/
│   ├── auth.py              # 認証モジュール
│   └── llm.py               # LLM ユーティリティ（Ollama / Groq）
└── data/                    # ローカルデータ置き場（.gitignore で除外）
```

---

## 🗺️ 今後の実装予定

- [ ] データ型自動判定（マイクロアレイ vs RNA-seq）
- [ ] 複数データセットのメタ解析
- [ ] 時系列データの対応
- [ ] レポート出力機能
- [ ] 研究室メンバー向けユーザー管理

---

## 📖 参考リソース

- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- [GEOparse ドキュメント](https://geoparse.readthedocs.io/)
- [gseapy ドキュメント](https://gseapy.readthedocs.io/)
- [STRING DB API](https://string-db.org/help/api/)
- [Streamlit ドキュメント](https://docs.streamlit.io/)
- [Ollama](https://ollama.com/)

---

## 👤 開発者

makoisaiah
