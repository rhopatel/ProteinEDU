import pygame
pygame.init()

screen = pygame.display.set_mode((500,500))
screen.fill((255,255,255))
pygame.display.set_caption('first game')
x = 50
y = 50
width = 40
height = 60
vel = 10
done = False

while not done:
    pygame.time.delay(100)
    screen.fill((255,255,255))

    for event in pygame.event.get():
        if (event.type == pygame.QUIT):
            done = True
            pygame.quit()
        if (event.type == pygame.MOUSEBUTTONDOWN):
            x,y = pygame.mouse.get_pos()
            print(x,y)
        
    pressed = pygame.key.get_pressed()
    if pressed[pygame.K_UP]: y -= vel
    if pressed[pygame.K_DOWN]: y += vel
    if pressed[pygame.K_LEFT]: x -= vel
    if pressed[pygame.K_RIGHT]: x += vel
    pygame.draw.rect(screen, (255,0,0), (x,y,width,height))
    pygame.display.update()

