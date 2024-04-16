SELECT ngs_forensic.marker_auto_strview.marker, 
	   ngs_forensic.marker_auto_strview.allele, 
       ngs_forensic.marker_auto_strview.sequence, 
       ngs_forensic.marker_auto_strview.count_seq, 
       ngs_forensic.marker_auto_strview.frequency,
       ngs_forensic.marker_auto_strview_flankingreg.sequence AS flank_sequence, 
       ngs_forensic.marker_auto_strview_flankingreg.count_seq AS flank_count_seq, ngs_forensic.marker_auto_strview_flankingreg.frequency AS flank_frequency
FROM ngs_forensic.marker_auto_strview
LEFT JOIN ngs_forensic.marker_auto_strview_flankingreg ON ngs_forensic.marker_auto_strview.allele = ngs_forensic.marker_auto_strview_flankingreg.allele 
AND ngs_forensic.marker_auto_strview.marker = ngs_forensic.marker_auto_strview_flankingreg.marker
AND ngs_forensic.marker_auto_strview_flankingreg.sequence LIKE CONCAT('%', ngs_forensic.marker_auto_strview.sequence, '%')
ORDER BY marker;