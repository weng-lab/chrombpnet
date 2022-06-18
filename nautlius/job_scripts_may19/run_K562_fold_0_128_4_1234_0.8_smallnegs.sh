apiVersion: batch/v1
kind: Job
metadata:
   name: chrombpnet.k562.128.4.1234.0.8.fold.0.smallnegs
   labels:
     jobgroup: jobexample
spec:
   template:
     metadata:
       name: jobexample
       labels:
         jobgroup: jobexample
     spec:
       affinity:
         nodeAffinity:
           preferredDuringSchedulingIgnoredDuringExecution:
           - weight: 1
             preference:
               matchExpressions:
               - key: gpu-type
                 operator: In
                 values:
                 - "titan-xp"
                 - "K40"
                 - "TITANRTX"
                 - "K40"
                 - "V100"
                 - "A100"
                 - "A40"
                 - "A40"
                 - "RTX6000"
                 - "3090"
       containers:
       - name: jobexample
         image: kundajelab/chrombpnet:dev 
         imagePullPolicy: Always
         resources:
           requests:
             cpu: 1
             memory: 16Gi
             nvidia.com/gpu: 1
           limits:
             cpu: 1
             memory: 32Gi
             nvidia.com/gpu: 1
         command:
         - /bin/bash
         - -c
         args:
         - nvidia-smi;
           mkdir /scratch/;
           cp -r /chrombpnet/chrombpnet/ /scratch/;
           cd /scratch/chrombpnet/;
           mkdir /scratch/chrombpnet/K562/; 
           cp -r /chrombpnet/ATAC_PE/K562/data/ /scratch/chrombpnet/K562/;
           cp -r /chrombpnet/reference/ /scratch/chrombpnet/;
           cp -r /chrombpnet/splits/ /scratch/chrombpnet/;
           bash atac_run_small_negs.sh 128 4 1234 0.8 K562 fold_0;
           cp -r /scratch/chrombpnet/models/* /chrombpnet/ATAC_PE/K562/; 
         volumeMounts:
         - mountPath: /chrombpnet
           name: chrombpnet
       restartPolicy: Never
       volumes:
         - name: chrombpnet
           persistentVolumeClaim:
             claimName: chrombpnet
       imagePullSecrets:
       - name: regcred

