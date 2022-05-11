import json
import telepot
import argparse
import time

# Create bot (see https://medium.com/@ManHay_Hong/how-to-create-a-telegram-bot-and-send-messages-with-python-4cf314d9fa3e):
# - On Telegram, search @ BotFather, send him a “/start” message
# - Send another “/newbot” message, then follow the instructions to setup a name and a username
# - Copy API token
# - Go to bot in Telegram, and press /start
# - Bot can be accessed by others as well.
# - Once activated, user_id will be added to getUpdates response

def send_bulletin(token,chat_id,bulletin,method):

	file = bulletin

	bot = telepot.Bot(token)

	# Check chat ID's
	# url = 'https://api.telegram.org/bot' + token + '/getUpdates'
	# resp = requests.get(url)
	# r_json = json.loads(resp.text)
	# print('r_json:')
	# print(r_json)

	# Method options are file and url
	if method == 'file':
		print(chat_id)
		bot.sendPhoto(chat_id, photo=open(file, 'rb'))
	elif method == 'document':
		print(chat_id)
		bot.sendDocument(chat_id, document=open(file, 'rb'))
	else:
		with open('bulletin.png', 'wb') as f:
			f.write(requests.get(file).content)
			f.close()
			time.sleep(3)

		print(chat_id)
		bot.sendPhoto(chat_id, photo=open('bulletin.png', 'rb'))

if __name__ == '__main__':

	# Get input from command line arguments
	parser = argparse.ArgumentParser(description = "Description for my parser")
	parser.add_argument("-T", "--token", help = "Telegram bot token", required = True, default = "")
	parser.add_argument("-C", "--chat_id", help = "Telegram chat ID", required = True, default = "")
	parser.add_argument("-B", "--bulletin", help = "Bulletin to be send", required = True, default = "")
	parser.add_argument("-M", "--method", help = "Specify file or URL as input", required = False, default = "url")

	argument = parser.parse_args()

	if argument.token:
		token = argument.token
		print('Bot token = ' + token)
	if argument.chat_id:
		chat_id = argument.chat_id
		print('Chat ID = ' + chat_id)
	if argument.bulletin:
		bulletin = argument.bulletin
		print('Bulletin filename = ' + bulletin)
	if argument.method:
		method = argument.method
		print('Method = ' + method)

	send_bulletin(token,chat_id,bulletin,method)
