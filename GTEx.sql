CREATE DATABASE `GTEX`;
SHOW DATABASES;

USE gtex;
create table first100gene(
gene_tissue_1 char(50),
gene_tissue_2 char(50),
bicor char(10),
pvalue char(20),
gene_symbol_1 char(30),
tissue_1 char(30),
gene_symbol_2 char(30),
tissue_2 char(30)
);


create table first1000gene(
gene_tissue_1 char(50),
gene_tissue_2 char(50),
bicor char(10),
pvalue char(20),
gene_symbol_1 char(30),
tissue_1 char(30),
gene_symbol_2 char(30),
tissue_2 char(30)
);
show tables;
select * from first100gene;
drop table first100gene;
select * from finalversion;

select * 
from finalversion
where tissue_1 = 'Adipose';

TRUNCATE TABLE finalversion;

FLUSH PRIVILEGES;

SHOW VARIABLES LIKE "secure_file_priv";
SHOW GLOBAL VARIABLES LIKE 'local_infile';
set global local_infile=true;


 
LOAD DATA local INFILE 'first 100 gene cocorrelation across two adipose tissues - melted data for SQL input.txt' INTO TABLE first100gene;
LOAD DATA INFILE 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\first 100 gene cocorrelation across two adipose tissues - melted data for SQL input.txt' INTO TABLE  first100gene;

BULK INSERT finalversion
FROM 'C:\Users\mingqiz7\Desktop\GTEx app\data\significant crosstissue enrichments ADIPOQ_Adipose - Subcutaneous - Final version.csv'
WITH
(
    FORMAT = 'csv',
    FIRSTROW = 2, -- as 1st one is header
    FIELDTERMINATOR = ',',  --CSV field delimiter
    ROWTERMINATOR = '\n',   --Use to shift the control to next row
    TABLOCK
)
GO

BULK INSERT finalversion
FROM 'C:\Users\mingqiz7\Desktop\GTEx app\data\significant crosstissue enrichments ADIPOQ_Adipose - Subcutaneous - Final version.csv'
WITH
(
        FORMAT='CSV',
        FIRSTROW=2
)
GRANT ALL PRIVILEGES ON *.* TO 'root'@'localhost' IDENTIFIED BY '1234567890zZ';
