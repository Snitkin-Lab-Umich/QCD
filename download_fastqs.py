import dropbox
import os

# Access Token
ACCESS_TOKEN = 'sl.CB4c1IgQgBz0rVS3fblM35GnZZkdQwi6pDTkF3h6VSkTc4YgaTADD2SvkC-MtF7sKBsO4jen62MytzQz4DaWZEWzZ4XIBSEkyWUwaOtbRmQKZhBjUjGW21Ro3MLyvKIs50qJioZLV-hJBhI'

# Create a Dropbox client
dbx = dropbox.Dropbox(ACCESS_TOKEN)

def download_files_from_dropbox(files_to_download):
    """Download multiple files from Dropbox to the local machine."""
    for dropbox_path, local_path in files_to_download:
        try:
            with open(local_path, "wb") as f:
                metadata, res = dbx.files_download(path=dropbox_path)
                f.write(res.content)
                print(f"Downloaded {dropbox_path} to {local_path}")
        except Exception as e:
            print(f"Error downloading {dropbox_path}: {e}")

# Define the file paths (Dropbox path and corresponding local path)
files_to_download = [
    ('/snitkin lab documents/test-pipelines/rush_kpc_110_r1.fastq.gz', '.test/dataRush_KPC_110_R1.fastq.gz'),
    ('/snitkin lab documents/test-pipelines/rush_kpc_110_r2.fastq.gz', '.test/data/Rush_KPC_110_R2.fastq.gz'),
    ('/snitkin lab documents/test-pipelines/rush_kpc_111_r1.fastq.gz', '.test/data/Rush_KPC_111_R1.fastq.gz'),
    ('/snitkin lab documents/test-pipelines/rush_kpc_111_r2.fastq.gz', '.test/data/Rush_KPC_111_R2.fastq.gz'),
    # Add more files as needed
]

# Download the files
download_files_from_dropbox(files_to_download)