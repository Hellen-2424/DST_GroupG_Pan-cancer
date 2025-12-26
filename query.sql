USE pancancer_immune;

SELECT
    pr.ligand,
    pr.receptor,
    pfL.category  AS ligand_category,
    pfR.category  AS receptor_category,
    COUNT(DISTINCT gc.cancer_type) AS num_cancers,
    AVG(gc.log2fc) AS mean_log2fc,
    MIN(gc.adj_p_value) AS best_p_value,
    GROUP_CONCAT(DISTINCT gc.cancer_type ORDER BY gc.cancer_type SEPARATOR ', ') AS cancer_types,
    GROUP_CONCAT(DISTINCT pw.pathway_name SEPARATOR '; ') AS associated_pathways
FROM pair AS pr
JOIN gene AS geneL
  ON geneL.gene_symbol = pr.ligand
LEFT JOIN protein_function AS pfL
  ON pfL.gene_id = geneL.gene_id
JOIN gene AS geneR
  ON geneR.gene_symbol = pr.receptor
LEFT JOIN protein_function AS pfR
  ON pfR.gene_id = geneR.gene_id
JOIN gene_cancer AS gc
  ON gc.gene_id IN (geneL.gene_id, geneR.gene_id)
LEFT JOIN pair_pathway AS pp
  ON pr.pair_id = pp.pair_id
LEFT JOIN pathway AS pw
  ON pp.pathway_id = pw.pathway_id
WHERE
    gc.cancer_type IN ('LUAD','COAD','BRCA')
    AND gc.adj_p_value < 0.10
GROUP BY
    pr.ligand,
    pr.receptor,
    pfL.category,
    pfR.category
HAVING
    COUNT(DISTINCT gc.cancer_type) >= 1
ORDER BY
    mean_log2fc DESC;
