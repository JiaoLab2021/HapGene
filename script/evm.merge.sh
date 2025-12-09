# Usage: bash evm.merge.sh $genome $gene_predictions $transcript_alignments $weight
if [ "$#" -ne 4 ]; then 
    echo "Usage: bash evm.merge.sh <genome> <gene_predictions> <protein_alignments> <weight>"
    exit 1
fi

genome=$1
gene_predictions=$2
protein_alignments=$3
weight=$4

EVM_BIN=$(command -v EVidenceModeler)

if [ -z "$EVM_BIN" ]; then
    echo "Error: EVidenceModeler not found in PATH. Please install and add to PATH."
    exit 1
fi

EVM_DIR=$(dirname "$EVM_BIN")
EvmUtils="$EVM_DIR/EvmUtils" 

perl ${EvmUtils}/partition_EVM_inputs.pl \
--genome $genome \
--segmentSize 1000000 \
--overlapSize 100000 \
--protein_alignments $protein_alignments \
--gene_predictions $gene_predictions \
--partition_listing partitions_list.out --partition_dir ./


perl ${EvmUtils}/write_EVM_commands.pl \
--genome $genome \
--gene_predictions $gene_predictions \
--protein_alignments $protein_alignments \
--weights $weight \
--partitions partitions_list.out \
--output_file_name evm.out \
> commands.list

parallel --jobs 30 < commands.list

perl ${EvmUtils}/recombine_EVM_partial_outputs.pl \
--partitions partitions_list.out \
--output_file_name evm.out

perl ${EvmUtils}/convert_EVM_outputs_to_GFF3.pl \
--partitions partitions_list.out \
--output evm.out \
--genome $genome

find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff3
