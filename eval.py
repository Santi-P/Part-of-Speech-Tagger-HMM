# simple evaluation script for HMM tagging
# Santi(chai) Pornavalai
# 1.4.19
# tested on python 3.7.2

import codecs

def evaluate(file1, file2):
    num = 0
    wrong = 0
    real =[]
    predicted = []
    tags = set()
    with codecs.open(file1,"r", encoding="iso8859-15") as f1, open(file2,"r",encoding="iso8859-15") as f2:
        
        for line1,line2 in zip(f1,f2):
            if len(line1)> 1 and len(line2) > 1:

                num += 1

                try:
                    word_1,tag_1 = line1.split()
                    word_2,tag_2 = line2.split()
                    real.append(tag_1)
                    tags.add(tag_1)
                    predicted.append(tag_2)
                    if tag_1 != tag_2:
                        wrong += 1
                except:
                    pass


                
                   # print(tag_1,tag_2)
    try:
        from sklearn.metrics import classification_report as report
        print("REPORT" ,report(real,predicted,labels=list(tags)))

    except ImportError:
        print("Sk-leanr module not found. skipping classification report")

    print( "accuracy",(num-wrong)/float(num), "%")



#evaluate("data/tiger_test.txt","results.txt")