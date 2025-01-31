#!/usr/bin/env python

import argparse
import json
import os
import shutil
import subprocess
import sys
from time import strftime, localtime

#
# Determine path to report template
#
top = os.getenv("KB_TOP")
report_deployed = os.path.join(top, "lib", "viral_assembly_report_template.html")
report_dev = os.path.join(top, "modules", "bvbrc_viral_assembly", "lib", "viral_assembly_report_template.html")
report_local = os.path.join("/home", "ac.mkuscuog", "git", "dev_container", "modules", "bvbrc_viral_assembly", "lib", "viral_assembly_report_template.html")
if os.path.exists(report_deployed):
  ASSEMBLY_REPORT_TEMPLATE = report_deployed
elif os.path.exists(report_dev):
  ASSEMBLY_REPORT_TEMPLATE = report_dev
else:
  ASSEMBLY_REPORT_TEMPLATE = report_local


DEFAULT_STRATEGY = "IRMA"
DEFAULT_IRMA_MODULE = "FLU"
MAX_RETRIES = 5
RETRY_DELAY = 10

def fetch_file_from_ws(ws_path, local_path):
  """Fetch a file from workspace to local path using p3-cp."""
  cmd = ["p3-cp", f"ws:{ws_path}", local_path]
  try:
    print(f"Fetching file from {ws_path} to {local_path}")
    subprocess.run(cmd, check=True)
    print(f"File successfully copied to {local_path}")
  except subprocess.CalledProcessError as e:
    print(f"Error fetching file from {ws_path}: {e}")
    return False
  return True

def fetch_fastqs_from_sra(sra_id, temp_dir="/tmp", output_dir="sra_fastqs"):
  os.makedirs(output_dir, exist_ok=True)
  cmd = ["fasterq-dump", "-t", temp_dir, "--outdir", output_dir, "--split-files", "-f", sra_id]

  for attempt in range(1, MAX_RETRIES + 1):
    print(f"Attempt {attempt}: Fetching FASTQs for SRA ID {sra_id} with command: {' '.join(cmd)}")
    try:
      subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
      print(f"Error fetching FASTQs for SRA ID {sra_id}: {e}")
      return None, None

    # Check if files exist after the command runs
    r1_path = os.path.join(output_dir, f"{sra_id}_1.fastq")
    r2_path = os.path.join(output_dir, f"{sra_id}_2.fastq")

    if os.path.exists(r1_path) or os.path.exists(r2_path):
      return (r1_path if os.path.exists(r1_path) else None,
              r2_path if os.path.exists(r2_path) else None)

    print(f"FASTQ files not found after attempt {attempt}. Retrying in {RETRY_DELAY} seconds...")
    time.sleep(retry_delay)

  print(f"No valid FASTQ files found for SRA ID {sra_id}.")
  return None, None

def concatenate_fasta_files(fasta_dir, output_fasta):
    """
    Concatenate all FASTA files in the output directory into one file.

    Args:
        fasta_dir (str): Directory containing output FASTA files.
        output_fasta (str): Path to the concatenated output FASTA file.
    """
    fasta_files = [
        os.path.join(fasta_dir, f)
        for f in os.listdir(fasta_dir)
        if f.endswith(".fasta") or f.endswith(".fa")
    ]

    if not fasta_files:
        print("No FASTA files found in the output directory.")
        return False

    print(f"Concatenating {len(fasta_files)} FASTA files into {output_fasta}...")
    with open(output_fasta, "w") as outfile:
        for fasta in fasta_files:
            with open(fasta, "r") as infile:
                outfile.write(infile.read())
    print("Concatenation complete.")
    return True

def run_quast(quast_output_dir, fasta_file, threads=12, min_contig=200):
    """
    Run QUAST on the concatenated FASTA file from IRMA output.

    Args:
        irma_output_dir (str): Directory containing IRMA output files.
        quast_output_dir (str): Directory to save QUAST output.
        threads (int): Number of threads to use for QUAST. Defaults to 12.
        min_contig (int): Minimum contig size for QUAST analysis. Defaults to 300.
    """
    os.makedirs(quast_output_dir, exist_ok=True)

    # Step 2: Run QUAST
    quast_cmd = [
        "quast.py",
        "-o", quast_output_dir,
        "-t", str(threads),
        "--min-contig", str(min_contig),
        fasta_file
    ]

    try:
        print(f"Running QUAST with command: {' '.join(quast_cmd)}")
        subprocess.run(quast_cmd, check=True)
        print("QUAST analysis completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running QUAST: {e}")

def run_irma(mode, input_file1, input_file2=None, output_dir="output"):
  irma_cmd = ["IRMA", mode, input_file1]
  if input_file2:
    irma_cmd.append(input_file2)
  irma_cmd.append(output_dir)

  try:
    print(f"Running IRMA with command: {' '.join(irma_cmd)}")
    subprocess.run(irma_cmd, check=True)
    print("IRMA run completed successfully!")
  except subprocess.CalledProcessError as e:
    print(f"Error running IRMA: {e}")

def get_software_version(software):
  try:
    result = subprocess.run(
      [software, "--version"],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True,
    )
    # Get only the first line of the output
    version = result.stdout.splitlines()[0]
    print(f"{software} version: {version}")
    return version.strip()  # Strip any extra spaces
  except Exception as e:
    print(f"Error fetching {software} version: {e}")
    return "Unknown"

def move_fasta_files(source_folder, destination_folder):
  for filename in os.listdir(source_folder):
    if filename.endswith(".fasta"):
      source_path = os.path.join(source_folder, filename)
      destination_path = os.path.join(destination_folder, filename)
            
      shutil.move(source_path, destination_path)
      print(f"Moved: {filename} -> {destination_folder}")

def generate_html_report(details):
    with open(ASSEMBLY_REPORT_TEMPLATE, "r") as template_file:
      template = template_file.read()

    # Set dynamic values
    report_date = strftime("%a, %d %b %Y %H:%M:%S", localtime())

    # Handle Quast section logic
    quast_html = details.get("quast_html", "")
    quast_txt_path = details.get("quast_txt", "")
    quast_section = ""
    if quast_html:
      quast_section = f"""
      <section>
          <h2>Quast Report</h2>
          <a href="{quast_html}">Quast HTML Report</a>
      """
      if quast_txt_path:
        with open(quast_txt_path, "r") as f:
          quast_txt_content = f.read()
        quast_section += f"""
        <pre class="preformatted">{quast_txt_content}</pre>
        """
      quast_section += "</section>"

    # Handle Tools table
    tools = details.get("tools", {})
    tools_table = "\n".join(
        f"<tr><td>{tool}</td><td>{version}</td></tr>"
        for tool, version in tools.items()
    )

    # Handle Errors
    errors_section = ""
    if "errors" in details and details["errors"]:
      '''
      errors_section = """
      <section>
          <h2>Errors</h2>
          <ul>
      """
      for error in details["errors"]:
        errors_section += f"<li>{error}</li>\n"
      errors_section += """
          </ul>
      </section>
      """
      '''
      errors_section = """
      <section style='color: #bb0505;'>
        <h2>Error</h2>
        <ul>
          <li>Assembly may have failed due to poor data quality or not having enough reads. We recommend running the FASTQC tool to assess read coverage and quality before retrying the assembly. You can access the FASTQC tool at <a href="https://www.bv-brc.org/app/FastqUtil" target="_blank">BV-BRC Fastq Utilities</a>, and the user guide is available <a href="https://www.bv-brc.org/docs/quick_references/services/fastq_utilities_service.html" target="_blank">here</a>.</li>
        </ul>
      </section>
      """

    # Replace placeholders
    html_content = template.replace("{{ report_date }}", report_date)
    html_content = html_content.replace("{{ quast_section }}", quast_section)
    html_content = html_content.replace("{{ tools_table }}", tools_table)
    html_content = html_content.replace("{{ errors_section }}", errors_section)

    return html_content

if __name__ == "__main__" :
  parser = argparse.ArgumentParser(description="Viral Assembly Script")
  parser.add_argument("-j", "--jfile", help="JSON file for job configuration", required=True)
  parser.add_argument("-o", "--output", help="Output directory. Defaults to current directory", required=False, default=".")

  args = parser.parse_args()

  #Load job data
  job_data = None
  try:
    with open(args.jfile, "r") as j:
      job_data = json.load(j)
  except Exception as e:
    print(f"Error loading JSON job file:\n {e}")
    traceback.print_exc(file=sys.stderr)
    sys.exit(-1)

  if not job_data:
    print("Job data is empty. Exiting.")
    sys.exit(-1) 

  print("Loaded job data:", job_data)  

  paired_end_lib = job_data.get("paired_end_lib", {})
  single_end_lib = job_data.get("single_end_lib", {})
  srr_id = job_data.get("srr_id", None)
  strategy = job_data.get("strategy", DEFAULT_STRATEGY)
  if strategy == "auto":
    strategy = DEFAULT_STRATEGY
  module = job_data.get("module", DEFAULT_IRMA_MODULE)

  # Ensure only one of the inputs is provided
  inputs_provided = sum(bool(x) for x in [paired_end_lib, single_end_lib, srr_id])
  if inputs_provided != 1:
    print("Error: Please provide exactly one of Paired End Library, Single End Librart, or SRR Id.")
    sys.exit(-1)

  # Setup output directory
  output_dir = args.output
  output_dir = os.path.abspath(output_dir)
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  os.chdir(output_dir)
  output_file = os.path.join(output_dir, job_data.get("output_file", "sample"))
  assembly_output_dir = os.path.join(output_dir, strategy)

  report_details = {
    "errors": []
  }

  # Process paired-end library
  if paired_end_lib:
    read1 = paired_end_lib.get("read1")
    read2 = paired_end_lib.get("read2")
    if not read1 or not read2:
      print("Error: Missing reads for paired-end library.")
      sys.exit(-1)

    # Fetch files locally
    local_read1 = os.path.join(output_dir, "read1.fasta")
    local_read2 = os.path.join(output_dir, "read2.fasta")
    if not (fetch_file_from_ws(read1, local_read1) and fetch_file_from_ws(read2, local_read2)):
      print("Error: Failed to fetch paired-end reads.")
      sys.exit(-1)

    run_irma(module, local_read1, local_read2, output_dir=assembly_output_dir)

  # Process single-end library
  elif single_end_lib:
    read = single_end_lib.get("read")
    if not read:
      print("Error: Missing read for single-end library.")
      sys.exit(-1)

    # Fetch file locally
    local_read = os.path.join(output_dir, "read.fasta")
    if not fetch_file_from_ws(read, local_read):
      print("Error: Failed to fetch single-end read.")
      sys.exit(-1)

    run_irma(module, local_read, output_dir=assembly_output_dir)

  # Process SRR ID
  elif srr_id:
    r1, r2 = fetch_fastqs_from_sra(srr_id, output_dir=output_dir)
    if not r1:
      print("Error: Failed to fetch FASTQs for SRR ID.")
      sys.exit(-1)

    run_irma(module, r1, r2, output_dir=assembly_output_dir)

  fasta_file = os.path.join(output_dir, output_file, f"{output_file}_all.fasta")
  is_concatenated = concatenate_fasta_files(assembly_output_dir, fasta_file)
  # Run QUAST on the concatenated FASTA file
  if is_concatenated:
    # Check if FASTA file is empty
    if os.path.getsize(fasta_file) == 0:
      print("FASTA file(s) generated by IRMA is empty.")
      report_details["errors"].append("FASTA file(s) generated by IRMA is empty.")
    else: 
      quast_dir = os.path.join(output_dir, "quast")
      try:
        run_quast(quast_dir, fasta_file)
        report_details["quast_txt"] = os.path.join(quast_dir, "report.txt")
        report_details["quast_html"] = "quast/report.html"
      except Exception:
        report_details["errors"].append("QUAST failed to run.")
  else:
     report_details["errors"].append("FASTA file was not generated by IRMA.")

  # Move individual segment files to the main output folder
  move_fasta_files(assembly_output_dir, output_dir) 

  # Retrieve software versions
  report_details["tools"] = {
    "IRMA": get_software_version("IRMA"),
    "quast": get_software_version("quast.py")
  }

  # Generate HTML report
  html_report = generate_html_report(report_details)
  with open(os.path.join(output_dir, "AssemblyReport.html"), "w") as output_file:
    output_file.write(html_report)
    
