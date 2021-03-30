

#input: integer money and array coins of d positive integers which are the denominations
#output: the min nimber of coins w demoninations coins that changes money

def dp_change(money, coins):
	MinNumCoins = [0]
	#print(MinNumCoins)

	for m in range(1, money + 1):
		MinNumCoins.append(1000) #said to go to infitiny, I arbitrarily picked this and it worked, but what is the real rule for how big I need to go?
		#print(MinNumCoins)
		for i in coins:
			if m >= i:
				if MinNumCoins[m-i] + 1 < MinNumCoins[m]:
					MinNumCoins[m] = MinNumCoins[m-i] + 1
	print(MinNumCoins)
	return MinNumCoins[money]



def main():
	amount = 40
	coins = [1, 5, 10, 20, 25, 50]

	min_num = str(dp_change(amount, coins))

	#Correct answer should be 2
	if min_num == '2':
		print('Correct you got 2')
	else: 
		print('Incorrect you got ' + min_num)


if __name__ == '__main__':
    main()