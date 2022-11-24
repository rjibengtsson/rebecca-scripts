import random

lst = []

for i in range(1,1000000):
    var1 = random.randint(1,4600000)
    var2 = random.randint(1,4600000)
    diff = abs(var1 - var2)
    lst.append(diff)


for i in range(1000,100000,1000):
    list1  = [k for k in lst if k < i]
    lst_len = len(list1)
    print(str(i), "=", lst_len)