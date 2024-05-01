import re

with open('periodicDipole.txt', 'r') as file:
    text = file.read()

text = re.sub(r'periodicDipolePoints=((?:.|\s)*?);', 'periodicDipolePoints=100;', text)
print(text)

with open('TEST.txt', 'w') as file:
    file.write(text)
