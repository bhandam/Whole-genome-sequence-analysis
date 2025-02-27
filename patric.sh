#take all the file from the spades to the patric folder
#to find teh genome that is close to
for spades in results/spades/AG*/contigs.fasta; do
sample_ID=$(dirname "$spades")
final_sample=$(basename "$sample_ID")
cp "$spades" results/patric/"$final_sample".fasta
done

#make directory for 12 serotypes
mkdir results/patric/{anatum,brenderup,cerro,dublin,enteritidis,\
give,litchfield,motevideo,muenchen,oranineburg,newport,paratyphi}

#Anatum:2
cp results/patric/AG21-005{3,4}.fasta results/patric/anatum

#Branderup:5
cp results/patric/AG21-006[1-5].fasta results/patric/brenderup

#cerro: 15
cp results/patric/AG21-0056.fasta results/patric/cerro
cp results/patric/AG21-0559.fasta results/patric/cerro
cp results/patric/AG21-056[0-3].fasta results/patric/cerro
cp results/patric/AG21-057[4-9].fasta results/patric/cerro
cp results/patric/AG21-058[0-2].fasta results/patric/cerro
ls results/patric/cerro/ | wc -l


#Dublin:1
cp results/patric/AG21-0051.fasta results/patric/dublin

#Enteritidis:2
cp results/patric/AG21-006[0,6].fasta results/patric/enteritidis

#Give:6
cp results/patric/AG21-056[4-8].fasta results/patric/give
cp results/patric/AG21-0572.fasta results/patric/give
rm 

#Litchfield:1
cp results/patric/AG21-0055.fasta results/patric/litchfield

#Montevideo:5
cp results/patric/AG21-005[0,7,8,9].fasta results/patric/motevideo
cp results/patric/AG21-0592.fasta results/patric/motevideo

#Muenchen:3
cp results/patric/AG21-05{89,90,91}.fasta results/patric/muenchen

#Newport:6
cp results/patric/AG21-058[3-8].fasta results/patric/newport

#Oranienburg:4
paratyphi
cp results/patric/AG21-05{69,70,71,73}.fasta results/patric/oranineburg

#paratyphi:1
cp results/patric/AG21-0052.fasta results/patric/paratyphi
