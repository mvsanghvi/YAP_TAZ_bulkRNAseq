import os
import boto3

AWS_ACCESS_KEY=os.getenv('AWS_ACCESS_KEY')
AWS_SECRET_KEY=os.getenv('AWS_SECRET_KEY')
ACCOUNT_ID=os.getenv('ACCOUNT_ID')

s3 = boto3.client(
    's3',
    aws_access_key_id=AWS_ACCESS_KEY,
    aws_secret_access_key=AWS_SECRET_KEY,
    aws_account_id=ACCOUNT_ID
)

s3_folders = ['WT_P62', 'Y2_P37_28', 'YT2_P37_7_17']
for folder in s3_folders:
    print(f'Downloading file: {folder}/{folder}_1.fq.gz')
    s3.download_file('yaptazko-bulk-rnaseq', f'{folder}/{folder}_1.fq.gz', f'./{folder}_1.fq.gz')
    print(f'Downloading file: {folder}/{folder}_2.fq.gz')
    s3.download_file('yaptazko-bulk-rnaseq', f'{folder}/{folder}_2.fq.gz', f'./{folder}_2.fq.gz')                        