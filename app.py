# OmicsCompass - メインアプリ
# Streamlit のマルチページアプリについて:
# https://docs.streamlit.io/develop/concepts/multipage-apps/overview

import streamlit as st
from components.auth import login

st.set_page_config(
    page_title="OmicsCompass",
    page_icon="🧭",
    layout="wide"
)

# ログイン確認
if not login():
    st.stop()

# ログイン後のメインページ
st.title("🧭 OmicsCompass")
st.subheader("オミクスデータから経路を探索するツール")

st.markdown("""
### ナビゲーション
左のサイドバーからページを選んでください。

| ページ | 内容 |
|--------|------|
| 📥 Data Fetch | GEO からRNA-seq データを取得 |
| 🔬 Analysis | 差分発現解析 |
| 📊 Visualization | パスウェイ解析と可視化 |
""")