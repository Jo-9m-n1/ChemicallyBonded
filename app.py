from flask import Flask, render_template, request, jsonify, session
from google import genai
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from dotenv import load_dotenv
import os
import re
import prototype_vsepr

load_dotenv()

WELCOME_MESSAGE = "Hi! Welcome to my office, how can I help you today?"

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


class Message(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    conversation_id = db.Column(db.Integer, db.ForeignKey("conversation.id"), nullable=False)
    is_user = db.Column(db.Boolean, nullable=False)
    text = db.Column(db.Text, nullable=False)
    lewis_structure = db.Column(db.String(500))
    vsepr = db.Column(db.String(500))

@app.route("/")
def index():
    session.pop('story', None)
    return render_template("index.html")

@app.route("/chat", methods=["POST"])
def chat():
    excuses = [
        "Dr. Lewis' wife just found a 'Tinder Platinum' charge on their joint card. He is currently sprinting toward the backyard. Try again later.",
        "He's been cornered near the shed. She has a heavy frying pan. Try again later.",
        "He's now hiding in the attic. He also forgot that it was their 10th anniversary. Try again later.",
        "The attic hatch has been breached. He is currently trying to negotiate a peace treaty. Try again later.",
        "He's been caught. He is currently begging for forgiveness in the driveway while neighbors watch. Try again later.",
        "His wife hit him with a frying pan after finding out that he texted someone named Michelle. Try again later.",
        "He's officially 'missing' after his wife went through his browser history. Try again later.",
        "You are invited to his funeral."
    ]

    data = request.json
    user_message = data.get("message")
    conv_id = data.get("conversation_id")
    verbosity = data.get("verbosity", "medium")

    conv = Conversation.query.filter_by(id=conv_id, user_id=1).first() if conv_id else None
    if not conv or not user_message:
        return jsonify({"reply": "Pick a chat from the list or start a new one."}), 400

    verbosity_map = {
        "not at all": "Provide only essential function calls and a brief greeting. Minimal text. If you are using LEWISFUNCTION or VSEPRFUNCTION, output ONLY the function call with no text at all.",
        "less": "Be extremely concise. Get straight to the point in 1-2 sentences.",
        "medium": "Provide a standard, professional explanation of moderate length.",
        "a lot": "Be very thorough. Explain the reasoning behind the concepts and provide detailed examples. Make sure that the student understand it."
    }

    selected_constraint = verbosity_map.get(verbosity, verbosity_map["medium"])

    LewisStructure = None
    VSEPRModel = None

    past_messages = Message.query.filter_by(conversation_id=conv_id).order_by(Message.id.asc()).all()

    history = []
    for msg in past_messages:
        role = "user" if msg.is_user else "model"
        history.append({"role": role, "parts": [{"text": msg.text}]})

    history.append({"role": "user", "parts": [{"text": user_message}]})

    try:
        with open("prompt.txt", "r", encoding="utf-8") as f:
            system_prompt = f.read()

        system_prompt_with_verbosity = f"{system_prompt}\n\nVERBOSITY GUIDELINE: {selected_constraint}"
        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=history,
            config={
                "system_instruction": system_prompt_with_verbosity
            }
        )

        reply = response.text

        CallLewis = re.search(r"LEWISFUNCTION\((.*?)\)", reply)
        Call3D = re.search(r"VSEPRFUNCTION\((.*?)\)", reply)

        if CallLewis:
            formula = CallLewis.group(1)
            response = prototype_vsepr.main(formula)
            if response is None:
                reply = re.sub(r"LEWISFUNCTION\(.*?\)", "", reply).strip()
                LewisStructure = f"static/generated/molecule_{formula}.svg"
            else:
                reply = response

        if Call3D:
            formula = Call3D.group(1)
            response = prototype_vsepr.main(formula)
            if response is None:
                reply = re.sub(r"VSEPRFUNCTION\(.*?\)", "", reply).strip()
                VSEPRModel = f"static/generated/vsepr_{formula}.html"
            else:
                reply = response

        db.session.add(Message(conversation_id=conv.id, is_user=True, text=user_message))
        db.session.add(
            Message(
                conversation_id=conv.id,
                is_user=False,
                text=reply,
                lewis_structure=LewisStructure,
                vsepr=VSEPRModel,
            )
        )
        db.session.commit()

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
    db.session.add(
        Message(conversation_id=new_chat.id, is_user=False, text=WELCOME_MESSAGE)
    )
    db.session.commit()
    session["conv_id"] = new_chat.id
    return jsonify({'id': new_chat.id})


@app.route("/conversation/<int:cid>/messages")
def conversation_messages(cid):
    conv = Conversation.query.filter_by(id=cid, user_id=1).first()
    if not conv:
        return jsonify([])
    rows = (
        Message.query.filter_by(conversation_id=cid).order_by(Message.id.asc()).all()
    )
    out = []
    for m in rows:
        out.append(
            {
                "role": "user" if m.is_user else "assistant",
                "text": m.text,
                "lewis_structure": m.lewis_structure,
                "vsepr": m.vsepr,
            }
        )
    return jsonify(out)

if __name__ == "__main__":
    with app.app_context():
        db.create_all()
    app.run(debug=True)