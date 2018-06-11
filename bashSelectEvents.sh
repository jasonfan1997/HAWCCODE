#!/bin/sh
#SBATCH --mem=500
#SBATCH -t 0-20:00:00

eval `$HAWCSOFT/setup.sh`
eval `/data/disk01/home/tcapistran/software/aerie/build/hawc-config --env-sh`
export CONFIG_HAWC=/data/disk01/home/tcapistran/software/config-hawc

DIRRUN="/data/disk01/home/tcapistran/Doctorado/HKU/Map/SelectEvts"
cd $DIRRUN
YEAR=2017
MONTHI=2
MONTHF=2

INDIR="/data/archive/hawcroot/data/hawc/reconstructed/hawcprod/v2.02.02/config-38125/reco_xcdf/"
OUTDIR="/data/scratch/userspace/tcapistran/RedNeuronal/RealData"

#Coordenadas del cangrejo
FUENTE="Crab"
AR=1.4591
DEC=0.384

#Coordenadas de Mrk 421 
#FUENTE="Mrk421"
#AR=2.8992
#DEC=0.6669

#Coordenadas de Mrk 451
#FUENTE="Mrk451"
#AR=3.5096
#DEC=0.6387

#Coordenadas de Mrk 501
#FUENTE="Mrk501"
#AR=4.4238
#DEC=0.6939
mkdir -p $OUTDIR
mkdir -p $OUTDIR/$FUENTE
#Solo dejar pasar un cuadrado alrededor del cualquier fuente (se empezo con la del cangrejo)
TOLE=0.12
for ((im=$MONTHI; im<=$MONTHF; im++))
do
	if(($im<10))
	then
		ARCH="0"$im
	else
		ARCH=$im
	fi
	if [ -d $INDIR/$YEAR/$ARCH ]
	then 
		ls $INDIR/$YEAR/$ARCH > lista_mes.txt
		for line_mes in $(cat lista_mes.txt)	
		do 
			mkdir -p $OUTDIR/$FUENTE/$YEAR/$ARCH	
			NOMSAL=$OUTDIR/$FUENTE/$YEAR/$ARCH"/all_rec_trig_"$line_mes".xcd"
			echo $NOMSAL
			if [ -f $NOMSAL ] || [ ! -d $INDIR/$YEAR/$ARCH/$line_mes ]
			then 
				echo "Ya existe el archivo o no existe el directorio "$NOMSAL
			else
				ls $INDIR/$YEAR/$ARCH/$line_mes/*.xcd > lista.txt
				ENTRADAS=""
				for line in $(cat lista.txt)
				do
						ENTRADAS=$ENTRADAS" "$line
				done
#				echo $ENTRADAS
#				xcdf-utility count  $ENTRADAS
#				echo $NOMSAL
				xcdf-utility select "(rec.coreFitStatus==0) && (rec.angleFitStatus == 0) && ((rec.ra-$AR) < $TOLE) && (($AR-rec.ra) < $TOLE) && ((rec.dec-$DEC) < $TOLE) && (($DEC-rec.dec) < $TOLE)" -o $NOMSAL $ENTRADAS

			fi
		done
	else 
		echo "No existe directorio "$INDIR/$YEAR/$ARCH
	fi
done
rm lista*
echo "Fin"
