#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --partition=jnovembre
#SBATCH --account=pi-jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=nb-tunnel
#SBATCH --output=nb-log-%J.out
#SBATCH --error=nb-log-%J.err

ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)
login=midway2.rcc.uchicago.edu


# These will be written to `nb-log-%J.out`
echo -e "\n\n   Copy & paste the following to your local computer's terminal to tunnel
to the remote server:\n"
echo        "   \`\`\`bash
"
echo        "   ssh -N -L $ipnport:$ipnip:$ipnport $USER@$login -v
"
echo        "   \`\`\`
"
echo -e "\n\n   Then use a web browser on your local computer to open the following
address:\n "
echo        "   \`\`\`google-chrome
"
echo        "   localhost:$ipnport
"
echo -e     "   \`\`\`\n\n"
echo        "   You may be asked to provide a token (a string of 48 letters + digits).
Please find it in the \`nb-log-*.err\` file under your current directory.
"
sleep 1
# Start jupyter
cd

module load python/3.7.0
module load plink
module load bcftools

jupyter-lab --no-browser --port=$ipnport --ip=$ipnip --log-level='ERROR'

# Then follow instructions above.
