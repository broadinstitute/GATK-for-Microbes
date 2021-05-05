FROM python:3

WORKDIR /usr/src/app

COPY ./requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY mapq0_cnv.py .

CMD [ "python", "./mapq0_cnv.py" ]