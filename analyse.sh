#!/bin/bash
if [ -L $0 ]; then
        BASE_DIR=`dirname $(readlink $0)`
else
        BASE_DIR=`dirname $0`
fi
basepath=$(cd $BASE_DIR; pwd)
JC_HOME=$basepath
TARGETS="/nfs/export3_25T/fep-rbfe/tyk2/parallel_test"

for target in $TARGETS;do
	cd $target
	for test in test1 test2 test3 test4 test5;do
		for pair in 23-55 43-55;do
			for system in complex ligands;do
                        	for edge in $(echo $system|cut -c 1)M2A $(echo $system|cut -c 1)M2B;do
                                	cd ${test}/${pair}/FEP/${system}/${edge}
					python /nfs/zli_gpu2/bin/developed_softwares/AlchemConvTools/one_end_fe_aly.py -i /nfs/zli_gpu2/bin/developed_softwares/AlchemConvTools/amber_aly.in
					cd $target
				done
			done
		done
	done
done
