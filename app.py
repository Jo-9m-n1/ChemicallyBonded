from flask import Flask, render_template, request, jsonify
from google import genai
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from dotenv import load_dotenv
import os
import re
import prototype_vsepr

load_dotenv()

app = Flask(__name__)

API_KEY = os.getenv('API')
client = genai.Client(api_key=API_KEY)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chats.db'
db = SQLAlchemy(app)

class Conversation(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, nullable=False)
    title = db.Column(db.String(200), nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.now)

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

        reply = response.text

        CallLewis = re.search(r"LEWISFUNCTION\((.*?)\)", reply)
        Call3D = re.search(r"VSEPRFUNCTION\((.*?)\)", reply)

        LewisStructure = None
        VSEPRModel = None


        if CallLewis:
            formula = CallLewis.group(1)
            
            prototype_vsepr.main(formula)
            
            reply = re.sub(r"LEWISFUNCTION\(.*?\)", "", reply).strip()
            

            LewisStructure = "/molecule.svg"

        elif Call3D:
            formula = Call3D.group(1)
            
            prototype_vsepr.main(formula)
            
            reply = re.sub(r"VSEPRFUNCTION\(.*?\)", "", reply).strip()
            

            VSEPRModel = "/static/html_visualisation.html" 

        return jsonify({
            "reply": reply,
            "LewisStructure": LewisStructure,
            "VSEPR": VSEPRModel,
        })
    
    except Exception as e:
        return jsonify({"reply": "Something went wrong."}), 500

@app.route('/get_conversations')
def get_conversations():
    temp_user_id = 1 
    convs = Conversation.query.filter_by(user_id=temp_user_id).order_by(Conversation.created_at.desc()).all()
    return jsonify([{'id': c.id, 'title': c.title, 'date': c.created_at.strftime('%Y-%m-%d')} for c in convs])

@app.route('/newchat', methods=["POST"])
def newchat():
    new_chat = Conversation(user_id = 1, title=f"{datetime.now().strftime('%H:%M:%S')}")
    db.session.add(new_chat)
    db.session.commit()
    return jsonify({'id': new_chat.id})

if __name__ == "__main__":
    with app.app_context():
        db.create_all()
    app.run(debug=True)