import os
import boto3

AWS_ACCESS_KEY=os.getenv('AWS_ACCESS_KEY')
AWS_SECRET_KEY=os.getenv('AWS_SECRET_KEY')
ACCOUNT_ID=os.getenv('ACCOUNT_ID')

s3= boto3.client(
    's3',
    aws_access_key_id=AWS_ACCESS_KEY,
    aws_secret_access_key=AWS_SECRET_KEY,
    aws_account_id=ACCOUNT_ID
)

s3.download_file('yaptazko-bulk-rnaseq', 'WT_P62/WT_P62_1.fq.gz', './WT_P62_1.fq.gz')