from muf import * 

try:
    obj=plot_IRI(version_IRI=2016)

    obj.run()
except:

    with open("/scripts/status.txt","w") as file:

        file.write("0")
else:

    with open("/scripts/status.txt","w") as file:

        file.write("1")
