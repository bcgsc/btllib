import btllib

with open("indexlr.fa", btllib.SeqReaderFlag.LONG_MODE) as reader:
	for record in reader:
		print(record.id)