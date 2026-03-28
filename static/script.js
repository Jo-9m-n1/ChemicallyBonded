const chat = document.getElementById('chat');
const input = document.getElementById('message-input');
const btn = document.getElementById('sendQuestion');

async function sendMessage() {
    const message = input.value.trim();
    if (!message) return;
    const userMsg = document.createElement('div');
    userMsg.className = 'messages messages--user';
    userMsg.textContent = message;
    chat.appendChild(userMsg);
    chat.scrollTop = chat.scrollHeight;
    input.value = '';
    const res = await fetch('/chat', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ message })
    });
    const data = await res.json();
    const botMsg = document.createElement('div');
    botMsg.className = 'messages messages--bot';
    botMsg.textContent = data.reply;
    chat.appendChild(botMsg);
    chat.scrollTop = chat.scrollHeight;
}

btn.addEventListener('click', sendMessage);
input.addEventListener('keypress', e => {
    if (e.key === 'Enter') sendMessage();
});
