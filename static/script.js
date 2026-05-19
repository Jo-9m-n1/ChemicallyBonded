const chat = document.getElementById('chat');
const input = document.getElementById('message-input');
const btn = document.getElementById('sendQuestion');
let currentConversationId = null;
const color1 = document.getElementById('color1');
const color2 = document.getElementById('color2');
const color3 = document.getElementById('color3');

async function sendMessage() {
    const message = input.value.trim();
    const verbosityDropdwon = document.getElementById('verbosity-select');
    const verbosity = verbosityDropdwon ? verbosityDropdwon.value: 'medium';

    if (!message) return;

    const userMsg = document.createElement('div');
    userMsg.className = 'messages messages--user';
    userMsg.textContent = message;
    chat.appendChild(userMsg);

    input.value = '';

    const typingBubble = document.createElement('div');
    typingBubble.className = 'messages messages--bot typing-indicator';
    typingBubble.innerHTML = '<span></span><span></span><span></span>'; 
    chat.appendChild(typingBubble);
    chat.scrollTop = chat.scrollHeight;

    try {
        const res = await fetch('/chat', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ 
                message, 
                conversation_id: currentConversationId, 
                verbosity: verbosity 
            })
        });

        const data = await res.json();

        if (chat.contains(typingBubble)) {
            chat.removeChild(typingBubble);
        }

        const timeDiv = document.createElement('div');
        timeDiv.style.cssText = "font-size: 0.67em; color: #747474; margin-left: 18px; margin-bottom: 1px; font-family: monospace;";
        timeDiv.textContent = `Processed in ${data.time || '0.000'}s`;

        chat.appendChild(timeDiv);

        if (data.reply && data.reply.trim().length > 0) {

            const botMsg = document.createElement('div');
            botMsg.className = 'messages messages--bot';
            botMsg.textContent = data.reply;
            chat.appendChild(botMsg);
        }

        if (data.LewisStructure || data.VSEPR) {
            const chemCard = createChemistryCard(data.VSEPR, data.LewisStructure);
            chat.appendChild(chemCard);
        }

        chat.scrollTop = chat.scrollHeight;

    } catch (e) {
        if (chat.contains(typingBubble)) {
            chat.removeChild(typingBubble);
        }
        console.error("Fetch error:", e);
    }
}

function createChemistryCard(vsepr, lewis) {
    const container = document.createElement('div');
    container.className = 'chemistry-card';
    
    container.innerHTML = `
        <div class="viewer-container">
            ${vsepr ? `<iframe src="/${vsepr}" class="mol-viewer"></iframe>` : ''}
            ${lewis ? `<div class="lewis-wrap"><img src="/${lewis}" alt="Lewis structure" class="lewis-structure"></div>` : ''}
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
            document.getElementById('message-input').value = '';
            await openConversation(data.id);
            loadHistory();
        }
    } catch (e) { console.error(e); }
}

async function loadHistory() {
    const res = await fetch('/get_conversations');
    const conversations = await res.json();
    const historyList = document.getElementById('history-list');

    historyList.innerHTML = conversations.map(c => `
        <div class="history-item" data-id="${c.id}">
            <strong>${c.title}</strong>
            <br><small>${c.date}</small>
        </div>
    `).join('');
}

async function openConversation(id) {
    const res = await fetch('/conversation/' + id + '/messages');
    const messages = await res.json();
    currentConversationId = Number(id);
    const chatEl = document.getElementById('chat');
    chatEl.innerHTML = '';
    if (!messages.length) {
        chatEl.innerHTML =
            '<div class="messages--bot messages">Hi! Welcome to my office, how can I help you today?</div>';
    } else {
        for (const m of messages) {
            if (m.role === 'user') {
                const div = document.createElement('div');
                div.className = 'messages messages--user';
                div.textContent = m.text;
                chatEl.appendChild(div);
            } else {
                const div = document.createElement('div');
                if (m.text && m.text.trim().length !== 0) { 
                    div.className = 'messages messages--bot';
                    div.textContent = m.text;
                    chatEl.appendChild(div);
                } if (m.lewis_structure || m.vsepr) {
                    chatEl.appendChild(createChemistryCard(m.vsepr, m.lewis_structure));
                }
            }
        }
    }
    chatEl.scrollTop = chatEl.scrollHeight;
}

document.addEventListener('DOMContentLoaded', () => {
    const historyList = document.getElementById('history-list');
    if (historyList) {
        historyList.addEventListener('click', (e) => {
            const item = e.target.closest('.history-item');
            if (!item || !item.dataset.id) return;
            openConversation(item.dataset.id);
        });
    }

    const homeView = document.getElementById('home-view');
    const profileView = document.getElementById('profile-view');
    const homeBtn = document.querySelector('.home'); 
    const profileBtn = document.querySelector('.profile');

    if (profileBtn) {
        profileBtn.addEventListener('click', () => {
            homeView.style.display = 'none';
            profileView.style.display = 'flex';
            homeBtn.classList.remove('tab-segment--active');
            profileBtn.classList.add('tab-segment--active');
        });
    }
    
    if (homeBtn) {
        homeBtn.addEventListener('click', () => {
            profileView.style.display = 'none';
            homeView.style.display = 'flex';
            profileBtn.classList.remove('tab-segment--active');
            homeBtn.classList.add('tab-segment--active');
        });
    }

    if (homeBtn) homeBtn.classList.add('tab-segment--active');
});

function updateBackground() {
    document.documentElement.style.setProperty('--bg-color-1', color1.value);
    document.documentElement.style.setProperty('--bg-color-2', color2.value);
    document.documentElement.style.setProperty('--bg-color-3', color3.value);
}

color1.addEventListener('input', updateBackground);
color2.addEventListener('input', updateBackground);
color3.addEventListener('input', updateBackground);

(async function init() {
    await startNewChat();
    await loadHistory();
})();