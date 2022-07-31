import pandas as pd
import seqlogo


def main():
    for i in range(1, 5):
        pwm = pd.read_table(f"PWM{i}.txt", delim_whitespace=True, header=None, comment='#')
        ppm = seqlogo.Ppm(pwm)
        seqlogo.seqlogo(ppm, ic_scale=False, format='svg', size='medium', filename=f'PWM{i}.svg')

if __name__ == '__main__':
    main()
