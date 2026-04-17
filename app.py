from flask import Flask, render_template, request, jsonify, session
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
app.secret_key = 'CanWeGet100Percent?'

class Conversation(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, nullable=False)
    title = db.Column(db.String(200), nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.now)

@app.route("/")
def index():
    session.pop('story', None)
    return render_template("index.html")

@app.route("/chat", methods=["POST"])
def chat():
    excuses = [
        "Dr. Lewis's wife just found a 'Tinder Platinum' charge on their joint card. He is currently sprinting toward the backyard.",
        "He's been cornered near the shed. She has a heavy frying pan.",
        "He's now hiding in the attic. He also forgot that it was their 10th anniversary.",
        "The attic hatch has been breached. He is currently trying to negotiate a peace treaty",
        "He's been caught. He is currently begging for forgiveness in the driveway while neighbors watch.",
        "His wife hit him with a frying pan after finding out that he texted someone named Michelle.",
        "He's officially 'missing' after his wife went through his browser history.",
        "You are invited to his funeral."
    ]

    data = request.json
    user_message = data.get("message")
    
    LewisStructure = None
    VSEPRModel = None

    try:
        with open("prompt.txt", "r", encoding="utf-8") as f:
            system_prompt = f.read()

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

        if CallLewis:
            formula = CallLewis.group(1)
            prototype_vsepr.main(formula)
            reply = re.sub(r"LEWISFUNCTION\(.*?\)", "", reply).strip()
            LewisStructure = "/molecule.svg"

        if Call3D:
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
        if 'story' not in session:
            session['story'] = 0
        
        current_index = session['story']
        
        if current_index >= len(excuses):
            current_index = len(excuses) - 1
            
        current_msg = excuses[current_index]
        session['story'] = min(current_index + 1, len(excuses) - 1)
        
        print(f"Error: {e}")
        return jsonify({"reply": current_msg}), 500

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