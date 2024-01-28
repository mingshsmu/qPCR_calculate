# qPCR_calculate
A repository containing some scripts for qPCR results calculating

Usage:
`Rscript.exe R_qPCR_analysis.R -i "2023-06-28 circTRHDE OE NOZ_20230629 131040.xlsx" -c NOZ_NC -r ACTB -o res`

Result:

**1. Expression:**
![Rplot_expression_NA](https://github.com/mingshsmu/qPCR_calculate/assets/47083119/c28fa64b-ff50-42eb-ab60-d24b8f16ce91)

**2. AmplicationCurve:**
![Rplot_AmplicationCurve_NA](https://github.com/mingshsmu/qPCR_calculate/assets/47083119/9d30789e-20db-4276-8efc-06e545e11063)

**3. MeltCurve:**
![Rplot_MeltCurve_NA](https://github.com/mingshsmu/qPCR_calculate/assets/47083119/9592a494-b421-46a2-91e7-692110a416b5)

**4. Values**
| **samples** | **gene**  | **ct**      | **ct_mean** | **ct_sd**   | **ct_iqr**  | **keep** | **dct**      | **ddct**     | **exp**     |
|-------------|-----------|-------------|-------------|-------------|-------------|----------|--------------|--------------|-------------|
| NOZ_NC      | ACTB      | 15.62303604 | 15.52503368 | 0.085421431 | 0.078335862 | TRUE     | 0.098002362  | 0.098002362  | 0.934325817 |
| NOZ_NC      | ACTB      | 15.46636432 | 15.52503368 | 0.085421431 | 0.078335862 | TRUE     | -0.058669361 | -0.058669361 | 1.041504708 |
| NOZ_NC      | ACTB      | 15.48570068 | 15.52503368 | 0.085421431 | 0.078335862 | TRUE     | -0.039333    | -0.039333    | 1.02763861  |
| NOZ_OE      | ACTB      | 15.27670725 | 15.25423384 | 0.03699369  | 0.032585142 | TRUE     | 0.022473405  | 0.022473405  | 0.984543322 |
| NOZ_OE      | ACTB      | 15.27445732 | 15.25423384 | 0.03699369  | 0.032585142 | TRUE     | 0.020223473  | 0.020223473  | 0.986079949 |
| NOZ_OE      | ACTB      | 15.21153696 | 15.25423384 | 0.03699369  | 0.032585142 | TRUE     | -0.042696878 | -0.042696878 | 1.030037512 |
| NOZ_NC      | circTRHDE | 27.82566848 | 27.96607725 | 0.123140482 | 0.11502269  | TRUE     | 12.3006348   | -0.14040877  | 1.102217371 |
| NOZ_NC      | circTRHDE | 28.05571386 | 27.96607725 | 0.123140482 | 0.11502269  | TRUE     | 12.53068018  | 0.08963661   | 0.939759428 |
| NOZ_NC      | circTRHDE | 28.01684941 | 27.96607725 | 0.123140482 | 0.11502269  | TRUE     | 12.49181573  | 0.050772159  | 0.965419479 |
