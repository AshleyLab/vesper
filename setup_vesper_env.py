import os
import subprocess
import sys

def run_command(command, check=True):
    """
    Runs a shell command and handles errors gracefully.
    """
    try:
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=check, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{command}' failed with exit code {e.returncode}")
        sys.exit(1)

def patch_cstag_lazy_import():
    """
    Patch cstag's to_pdf.py to lazy-import weasyprint only when needed.
    """
    to_pdf_path = os.path.join(sys.prefix, "lib", "python3.9", "site-packages", "cstag", "to_pdf.py")

    if not os.path.exists(to_pdf_path):
        print(f"⚠️  cstag's to_pdf.py not found at {to_pdf_path}. Skipping patch.")
        return

    with open(to_pdf_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    import_removed = False
    for line in lines:
        if "from weasyprint import HTML" in line and not import_removed:
            print("🔧 Removing top-level weasyprint import...")
            import_removed = True
            continue
        new_lines.append(line)

    for i, line in enumerate(new_lines):
        if line.strip().startswith("def to_pdf("):
            indent = ' ' * (len(line) - len(line.lstrip()) + 4)
            new_lines.insert(i + 1, f"{indent}from weasyprint import HTML  # lazy import\n")
            print("✅ Injected lazy import inside `to_pdf()`.")
            break
    else:
        print("❌ Couldn't find `def to_pdf()` in cstag. No changes made.")
        return

    with open(to_pdf_path, "w") as f:
        f.writelines(new_lines)

    print("🎉 Successfully patched cstag to lazy-load weasyprint.")

def main():
    ### ensure script is running in a virtual environment ###
    if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
        print("Please activate your virtual environment before running this script.")
        sys.exit(1)

    ### Define installation path ###
    env_bin = os.path.join(sys.prefix, "bin")

    print(f"Virtual environment detected: {sys.prefix}")
    print(f"Installing tools into: {env_bin}")

    ### Upgrade pip and install Python dependencies ###
    run_command("pip install --upgrade pip")
    run_command("pip install -r requirements.txt")

    ### Install Samtools ###
    run_command("wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2")
    run_command("tar -xvjf samtools-1.16.1.tar.bz2")
    run_command("cd samtools-1.16.1 && ./configure --prefix=$VIRTUAL_ENV --disable-lzma --without-curses && make && make install && cd ..")

    ### Patch cstag's PDF import ###
    patch_cstag_lazy_import()

    ### Clean up temporary files ###
    run_command("rm -f samtools-1.16.1.tar.bz2")
    run_command("rm -rf samtools-1.16.1")

    ### Verify installations ###
    run_command(f"{env_bin}/samtools --version")

    print("\nAll tools installed successfully!")

if __name__ == "__main__":
    main()