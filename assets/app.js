const cookieStorage = {
    getItem: (key) => {
	const cookies = document.cookie
	      .split(';')
	      .map(cookie => cookie.split('='))
	      .reduce((acc, [key, value]) => ({ ...acc, [key.trim()]: value }), {});
	return cookies[key];
    },
    setItem: (key, value) => {
	document.cookie = `${key}=${value}`
    },
};

const storageType = cookieStorage;
const consentPropertyName = 'pcr_consent';

const shouldShowPopup = () => !storageType.getItem(consentPropertyName);
const saveToStorage = () => storageType.setItem(consentPropertyName, true);

window.onload = () => {
    const consentPopup = document.getElementById('consent-popup');
    const acceptBtn = document.getElementById('accept');

    const acceptFn = event => {
        event.preventDefault();
	saveToStorage(storageType);
	consentPopup.classList.add('hidden');
    };

    acceptBtn.addEventListener('click', acceptFn);
    
    if (shouldShowPopup(storageType)) {
	setTimeout(() => {
            consentPopup.classList.remove('hidden');
	}, 2000);
    }
};
