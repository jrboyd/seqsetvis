#!/bin/bash
cat /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/broad/H7_H3K27ME3_pooled_peaks.broadPeak | awk '{if ($1 == "chr19") print $0}' > H7_H3K27ME3_pooled_peaks.chr19.broadPeak
cat /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K27ME3_pooled/broadpeak/MCF7_H3K27ME3_pooled_peaks.broadPeak | awk '{if ($1 == "chr19") print $0}' > MCF7_H3K27ME3_pooled_peaks.chr19.broadPeak
cat /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K4ME3_pooled/MCF7_H3K4ME3_pooled_peaks.narrowPeak | awk '{if ($1 == "chr19") print $0}' > MCF7_H3K4ME3_pooled_peaks.chr19.narrowPeak
cat /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_peaks.narrowPeak | awk '{if ($1 == "chr19") print $0}' > H7_H3K4ME3_pooled_peaks.chr19.narrowPeak
cat /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_peaks.narrowPeak | awk '{if ($1 == "chr19") print $0}' > H7_H3K27ME3_pooled_peaks.chr19.narrowPeak
cat /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K27ME3_pooled/MCF7_H3K27ME3_pooled_peaks.narrowPeak | awk '{if ($1 == "chr19") print $0}' > MCF7_H3K27ME3_pooled_peaks.chr19.narrowPeak

ln -s /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/broad/H7_H3K27ME3_pooled_peaks.broadPeak full/
ln -s /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K27ME3_pooled/broadpeak/MCF7_H3K27ME3_pooled_peaks.broadPeak full/
ln -s /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K4ME3_pooled/MCF7_H3K4ME3_pooled_peaks.narrowPeak full/
ln -s /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_peaks.narrowPeak full/
ln -s /slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_peaks.narrowPeak full/
ln -s /slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K27ME3_pooled/MCF7_H3K27ME3_pooled_peaks.narrowPeak full/
