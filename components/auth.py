# OmicsCompass - 認証モジュール
# Streamlit Secrets を使った認証について:
# https://docs.streamlit.io/develop/concepts/connections/secrets-management
# streamlit-authenticator 0.4.x の使い方:
# https://github.com/mkhorasani/Streamlit-Authenticator

import streamlit as st
import streamlit_authenticator as stauth
import yaml

def load_auth():
    """Streamlit Secrets から認証情報を読み込む"""
    # ローカル開発時は auth_config.yaml を使い、
    # デプロイ時は Streamlit Cloud の Secrets を使う
    try:
        # Streamlit Cloud の Secrets から読み込む
        config = {
            "credentials": dict(st.secrets["credentials"]),
            "cookie": dict(st.secrets["cookie"])
        }
        # usernames の中身も dict に変換
        config["credentials"]["usernames"] = {
            k: dict(v) for k, v in st.secrets["credentials"]["usernames"].items()
        }
    except Exception:
        # ローカル開発時は yaml ファイルから読み込む
        with open("auth_config.yaml") as file:
            config = yaml.safe_load(file)

    authenticator = stauth.Authenticate(
        config["credentials"],
        config["cookie"]["name"],
        config["cookie"]["key"],
        config["cookie"]["expiry_days"],
    )
    return authenticator

def login():
    """ログイン画面を表示して認証状態を返す"""
    authenticator = load_auth()

    authenticator.login()

    if st.session_state.get("authentication_status") is False:
        st.error("ユーザー名またはパスワードが違います")
        return False
    elif st.session_state.get("authentication_status") is None:
        st.warning("ユーザー名とパスワードを入力してください")
        return False
    elif st.session_state.get("authentication_status"):
        authenticator.logout(location="sidebar")
        st.sidebar.success(f'ログイン中: {st.session_state["name"]}')
        return True