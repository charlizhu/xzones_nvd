def test():

    file = open("connect2.txt",'r')

    previousitems = []
    currentitems = []
    trackeditems = []
    maxcount = 0 # may or may not serve an actual purpose...

    for eachindex, eachrow in enumerate(file):

        eachitem = eachrow.split()

        # print(eachitem)

        if eachindex == 0:
            previousitems = eachitem
            pass

        currentitems = eachitem

        for thisitem in currentitems:
            if thisitem in previousitems:
                if thisitem not in trackeditems:
                    trackeditems.append(thisitem)
        for thatitem in trackeditems:
            if thatitem not in currentitems:
                trackeditems.pop(trackeditems.index(thatitem))
        # print(currentitems,previousitems)
        if len(set(currentitems).intersection(set(previousitems))) > 0:
            maxcount = maxcount + 1
        else:
            maxcount = 0

        previousitems = currentitems
        print(trackeditems, maxcount)

test()


