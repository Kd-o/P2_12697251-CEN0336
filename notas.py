#!/usr/bin/env python3

TOTAL = 0
CONTADOR_NOTAS = 0

while CONTADOR_NOTAS <= 10:
    k = int(input('Insira uma nota:'))
    TOTAL += k
    CONTADOR_NOTAS += 1

media = TOTAL/10
print('A media eh:',media)

import hashlib

with open("notas.py", "r") as f:
  script = f.read()

md5sum = hashlib.md5(script.encode("utf-8")).hexdigest()

print(md5sum)

