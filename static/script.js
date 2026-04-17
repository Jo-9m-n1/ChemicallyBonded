const chat = document.getElementById('chat');
const input = document.getElementById('message-input');
const btn = document.getElementById('sendQuestion');
const color1 = document.getElementById('color1');
const color2 = document.getElementById('color2');
const color3 = document.getElementById('color3');

async function sendMessage() {
    const message = input.value.trim();

    if (!message) return;

    const userMsg = document.createElement('div');
    userMsg.className = 'messages messages--user';
    userMsg.textContent = message;
    chat.appendChild(userMsg);
    chat.scrollTop = chat.scrollHeight;
    input.value = '';

    try {
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

        if (data.LewisStructure || data.VSEPR) {
            const chemCard = createChemistryCard(data.VSEPR, data.LewisStructure);
            chat.appendChild(chemCard);
        }

        chat.scrollTop = chat.scrollHeight;

    } catch (e) {
        console.error("Fetch error:", e);
    }
}

function createChemistryCard(vsepr, lewis) {
    const container = document.createElement('div');
    container.className = 'chemistry-card';
    
    container.innerHTML = `
        <div class="viewer-container">
            ${vsepr ? `<iframe src="/static/html_visualisation.html" class="mol-viewer"></iframe>` : ''}
            ${lewis ? `<object data="/static/molecule.svg" type="image/svg+xml" class="lewis-structure"></object>` : ''}
        </div>
    `;
    
    return container;
}

btn.addEventListener('click', sendMessage);
input.addEventListener('keypress', e => {
    if (e.key === 'Enter') sendMessage();
});

async function startNewChat() {
    try {
        const res = await fetch('/newchat', { method: 'POST' });
        const data = await res.json();

        if (data.id) {
            document.getElementById('chat').innerHTML = '<div class="messages--bot messages">Hi! Welcome to my office, how can I help you today?</div>';
            document.getElementById('message-input').value = '';
            loadHistory(); 
        }
    } catch (e) { console.error(e); }
}

async function loadHistory() {
    const res = await fetch('/get_conversations');
    const conversations = await res.json();
    const historyList = document.getElementById('history-list');

    historyList.innerHTML = conversations.map(chat => `
        <div class="history-item">
            <strong>${chat.title}</strong>
            <br><small>${chat.date}</small>
        </div>
    `).join('');
}

document.addEventListener('DOMContentLoaded', () => {
    const homeView = document.getElementById('home-view');
    const profileView = document.getElementById('profile-view');
    const homeBtn = document.querySelector('.home'); 
    const profileBtn = document.querySelector('.profile');

    if (profileBtn) {
        profileBtn.addEventListener('click', () => {
            homeView.style.display = 'none';
            profileView.style.display = 'flex';
        });
    }
    
    if (homeBtn) {
        homeBtn.addEventListener('click', () => {
            profileView.style.display = 'none';
            homeView.style.display = 'flex';
        });
    }
});

function updateBackground() {
    document.documentElement.style.setProperty('--bg-color-1', color1.value);
    document.documentElement.style.setProperty('--bg-color-2', color2.value);
    document.documentElement.style.setProperty('--bg-color-3', color3.value);
}

color1.addEventListener('input', updateBackground);
color2.addEventListener('input', updateBackground);
color3.addEventListener('input', updateBackground);

startNewChat()
window.onload = loadHistory;