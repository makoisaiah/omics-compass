# OmicsCompass - LLM ユーティリティ
# Ollama の Python クライアント:
# https://github.com/ollama/ollama-python
# Groq API:
# https://console.groq.com/docs/quickstart

import requests
import json

OLLAMA_URL = "http://localhost:11434"
OLLAMA_MODEL = "mistral"

def is_ollama_available():
    """Ollama がローカルで起動しているか確認"""
    try:
        response = requests.get(f"{OLLAMA_URL}/api/tags", timeout=2)
        return response.status_code == 200
    except Exception:
        return False

def ask_ollama(prompt: str) -> str:
    """ローカルの Ollama に問い合わせる"""
    response = requests.post(
        f"{OLLAMA_URL}/api/generate",
        json={
            "model": OLLAMA_MODEL,
            "prompt": prompt,
            "stream": False
        },
        timeout=60
    )
    return response.json()["response"]

def ask_groq(prompt: str, api_key: str) -> str:
    """Groq API に問い合わせる（フォールバック）"""
    response = requests.post(
        "https://api.groq.com/openai/v1/chat/completions",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        },
        json={
            "model": "llama-3.1-8b-instant",
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 500
        },
        timeout=30
    )
    return response.json()["choices"][0]["message"]["content"]

def ask_llm(prompt: str, groq_api_key: str = None) -> tuple[str, str]:
    """
    LLM に問い合わせる。ローカル優先、なければ Groq にフォールバック。
    戻り値: (回答テキスト, 使用したバックエンド名)
    """
    if is_ollama_available():
        return ask_ollama(prompt), "ローカル Mistral"
    elif groq_api_key:
        return ask_groq(prompt, groq_api_key), "Groq API"
    else:
        raise Exception("Ollama も Groq API キーも利用できません")
