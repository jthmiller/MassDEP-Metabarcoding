ref-dbs.sh


cp ~/old-home/watts/ref-database/rbcl/diat_barcode_v10_263bp-sklearn-classifier.qza ../../refdbs/rbcl/
cp ~/old-home/watts/ref-database/rbcl/diat_barcode_v10_263bp-seqs.qza ../../refdbs/rbcl/
cp ~/old-home/watts/ref-database/rbcl/diat_barcode_v10_263bp-tax.qza ../../refdbs/rbcl/


# train without the weights
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ../../refdbs/rbcl/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy ../../refdbs/rbcl/diat_barcode_v10_263bp-tax.qza \
  --o-classifier ../../refdbs/rbcl/diat_barcode_v10_263bp-sklearn-classifier_1.4.2.qza
