int n=0;

void setup() 
{
  Serial.begin(9600);
  pinMode(10, INPUT);
  pinMode(11, INPUT);


}


void loop() {
  n+=1;
  if((digitalRead(10)==1) || (digitalRead(11)==1)){
    Serial.println("!");
  }
  else{
    Serial.println(analogRead(A0));
  }
  delay(1);

  if (n==2850){
    Serial.end();
 
  }
}
