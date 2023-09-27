set search_path to soilweb, public;
\timing

-- 500ms
CREATE TEMP TABLE mukey_area AS
SELECT mukey, SUM(area_ac) AS area_ac
FROM soilweb.mukey_mlra_overlap
GROUP by mukey;


CREATE INDEX ON mukey_area (mukey);
VACUUM ANALYZE mukey_area;

-- 5.6 seconds
CREATE TEMP TABLE stock_fractions AS
SELECT stock_component_data.mukey, 
comppct_r * area_ac::numeric as area_fraction, 
oc_kg_sq_m, oc_kg_sq_m030, oc_kg_sq_m0100 
FROM stock_component_data 
JOIN mukey_area USING (mukey)
WHERE area_ac IS NOT NULL and area_ac > 0 ;

CREATE INDEX ON mukey_area (mukey);
VACUUM ANALYZE mukey_area;


-- convert ac -> m^2 4046.86
-- convert kg -> metric ton 0.001
SELECT 
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m * 0.001)) as soc_mton,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m030 * 0.001)) as soc_mton_030,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m0100 * 0.001)) as soc_mton_0100
FROM stock_fractions
WHERE mukey = '623418';

-- convert ac -> m^2 4046.86
-- convert kg -> Tg 1e-9
-- 3 seconds
SELECT 
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m * 1e-9)) as soc_Tg,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m030 * 1e-9)) as soc_Tg_030,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m0100 * 1e-9)) as soc_Tg_0100
FROM stock_fractions;

 soc_tg | soc_tg_030 | soc_tg_0100
--------+------------+-------------
  82625 |      34402 |       64744


-- convert ac -> m^2 4046.86
-- convert kg -> Gt 1/1000000000000
-- 3 seconds
SELECT 
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m * 1e-12), 2) as soc_Gt,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m030 * 1e-12), 2) as soc_Gt_030,
ROUND(SUM(area_fraction * 4046.86 * oc_kg_sq_m0100 * 1e-12), 2) as soc_Gt_0100
FROM stock_fractions;

 soc_gt | soc_gt_030 | soc_gt_0100
--------+------------+-------------
  82.63 |      34.40 |       64.74





