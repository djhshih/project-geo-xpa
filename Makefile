
all: GSE100855/raw GSE55484/raw GSE70818/raw
	

GSE100855/raw:
	cd GSE100855 && get.sh

GSE55484/raw:
	cd GSE55484 && get.sh

GSE70818/raw:
	cd GSE70818 && get.sh

