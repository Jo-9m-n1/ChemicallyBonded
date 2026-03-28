from flask import Flask, render_template, request, jsonify
from google import genai
from dotenv import load_dotenv
import os

load_dotenv()

app = Flask(__name__)

API_KEY = os.getenv('API')
client = genai.Client(api_key=API_KEY)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/chat", methods=["POST"])
def chat():
    data = request.json
    user_message = data.get("message")
    
    with open("prompt.txt", "r", encoding="utf-8") as f:
        system_prompt = f.read()

    try:
        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=user_message,
            config={
                "system_instruction": system_prompt
            }
        )
        return jsonify({"reply": response.text})
    except Exception as e:
        return jsonify({"reply": "Something went wrong."}), 500

if __name__ == "__main__":
    app.run(debug=True)