import os
import petl as etl
import pandas as pd
from sqlalchemy import create_engine, text

# Initialize database with performance optimizations
engine = create_engine('sqlite:///quantms.db')

# Configure SQLite for better performance
with engine.begin() as conn:
    conn.execute(text("PRAGMA journal_mode = WAL"))
    conn.execute(text("PRAGMA synchronous = NORMAL"))
    conn.execute(text("PRAGMA cache_size = -10000"))  # 10MB cache

# Dictionary to map file suffixes to their types
file_types = {
    'ms_info.tsv': 'ms_info',
    'msgf.tsv': 'msgf',
    'msgf_feat_perc.tsv': 'msgf_feat_perc',
    'msgf_feat_perc_pep.tsv': 'msgf_feat_perc_pep',
    'msgf_feat_perc_pep_filter.tsv': 'msgf_feat_perc_pep_filter',
    'msgf_feat_perc_pep_filter_psm.tsv': 'msgf_feat_perc_pep_filter_psm',
    'spectrum_df.tsv': 'spectrum'
}

# Dictionary to store all DataFrames by file type and sample
data_frames = {file_type: {} for file_type in file_types.values()}

# Process each file
for filename in os.listdir("./"):
    if not filename.endswith('.tsv'):
        continue

    # Find which file type this is
    file_type = None
    for suffix, f_type in file_types.items():
        if filename.endswith(suffix):
            file_type = f_type
            break

    if not file_type:
        print(f"Skipping {filename} - unrecognized file type")
        continue

    # Extract sample name (format: RA_10_1, RA_11_3, etc.)
    sample_name = '_'.join(filename.split('_')[:3])

    print(f"Processing {file_type} for sample {sample_name} ({filename})")

    try:
        # Read the file
        etl_table = etl.fromtsv(filename)

        # Add sample name as a column
        etl_table = etl.addfield(etl_table, 'sample_name', sample_name)

        # Convert to DataFrame
        df = etl.todataframe(etl_table)

        # Remove duplicates based on all columns within each sample's DataFrame
        df_no_duplicates = df.drop_duplicates()

        # Store in our dictionary
        data_frames[file_type][sample_name] = df_no_duplicates

        print(f"Successfully processed {filename} (duplicates removed)")

    except Exception as e:
        print(f"Error processing {filename}: {e}")

# Create tables in the database for each file type
for file_type, samples in data_frames.items():
    if not samples:
        print(f"No data found for {file_type}")
        continue

    # Combine all samples for this file type
    combined_df = pd.concat(samples.values(), ignore_index=True)

    # Write to SQL database
    try:
        with engine.begin() as conn:
            combined_df.to_sql(file_type, conn, if_exists='replace', index=False)
            print(f"Created/updated {file_type} table with {len(combined_df)} rows")
    except Exception as e:
        print(f"Error writing {file_type} table to database: {e}")

# Rename columns and create combined_score table using SQL
try:
    with engine.begin() as conn:
        # Rename columns using SQL operations
        print("Renaming columns in database tables...")

        # msgf_feat_perc_pep_filter: score → qvalue_score
        if 'msgf_feat_perc_pep_filter' in data_frames:
            conn.execute(text("""
                CREATE TABLE msgf_feat_perc_pep_filter_renamed AS
                SELECT *, score AS qvalue_score FROM msgf_feat_perc_pep_filter
            """))
            conn.execute(text("DROP TABLE msgf_feat_perc_pep_filter"))
            conn.execute(text("ALTER TABLE msgf_feat_perc_pep_filter_renamed RENAME TO msgf_feat_perc_pep_filter"))
            print("Renamed 'score' to 'qvalue_score' in msgf_feat_perc_pep_filter")

        # msgf: score → msgfplus_score
        if 'msgf' in data_frames:
            conn.execute(text("""
                CREATE TABLE msgf_renamed AS
                SELECT *, score AS msgfplus_score FROM msgf
            """))
            conn.execute(text("DROP TABLE msgf"))
            conn.execute(text("ALTER TABLE msgf_renamed RENAME TO msgf"))
            print("Renamed 'score' to 'msgfplus_score' in msgf")

        # msgf_feat_perc: score → percolator_score
        if 'msgf_feat_perc' in data_frames:
            conn.execute(text("""
                CREATE TABLE msgf_feat_perc_renamed AS
                SELECT *, score AS percolator_score FROM msgf_feat_perc
            """))
            conn.execute(text("DROP TABLE msgf_feat_perc"))
            conn.execute(text("ALTER TABLE msgf_feat_perc_renamed RENAME TO msgf_feat_perc"))
            print("Renamed 'score' to 'percolator_score' in msgf_feat_perc")

        # Create indexes for faster joining
        print("Creating indexes on join columns...")
        conn.execute(text("CREATE INDEX IF NOT EXISTS idx_msgf_join ON msgf(sample_name, mz, charge, sequence, start, end)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS idx_perc_join ON msgf_feat_perc(sample_name, mz, charge, sequence, start, end)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS idx_qvalue_join ON msgf_feat_perc_pep_filter(sample_name, mz, charge, sequence, start, end)"))

        # Perform the join entirely in SQL with the specified columns and remove duplicates
        print("Creating combined_score table via SQL join (removing duplicates)...")
        conn.execute(text("DROP TABLE IF EXISTS combined_score"))
        conn.execute(text("""
            CREATE TABLE combined_score AS
            SELECT DISTINCT
                m.sample_name,
                m.rt,
                m.mz,
                m.charge,
                m.aa_before,
                m.aa_after,
                m.sequence,
                m.start,
                m.end,
                m.protein_references,
                m.accessions,
                m.msgfplus_score,
                p.percolator_score,
                q.qvalue_score
            FROM msgf m
            INNER JOIN msgf_feat_perc p
                ON m.sample_name = p.sample_name
                AND m.mz = p.mz
                AND m.charge = p.charge
                AND m.sequence = p.sequence
                AND m.start = p.start
                AND m.end = p.end
                AND m.aa_before = p.aa_before
                AND m.aa_after = p.aa_after
                AND m.protein_references = p.protein_references
                AND m.accessions = p.accessions
                AND ABS(m.rt - p.rt) < 0.1 -- Adjust tolerance as needed
            INNER JOIN msgf_feat_perc_pep_filter q
                ON m.sample_name = q.sample_name
                AND m.mz = q.mz
                AND m.charge = q.charge
                AND m.sequence = q.sequence
                AND m.start = q.start
                AND m.end = q.end
                AND m.aa_before = q.aa_before
                AND m.aa_after = q.aa_after
                AND m.protein_references = q.protein_references
                AND m.accessions = q.accessions
                AND ABS(m.rt - q.rt) < 0.1 -- Adjust tolerance as needed
        """))

        # Verify the result
        count = conn.execute(text("SELECT COUNT(*) FROM combined_score")).scalar()
        print(f"Successfully created combined_score table with {count} rows (duplicates removed)")
        print("Columns in combined_score table:")
        cols = conn.execute(text("PRAGMA table_info(combined_score)")).fetchall()
        for col in cols:
            print(f"- {col[1]}")  # Column name is at index 1

except Exception as e:
    print(f"Error during SQL operations: {e}")
finally:
    engine.dispose()

print("Processing complete!")