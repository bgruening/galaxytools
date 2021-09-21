import numpy as np
import tensorflow as tf


(mnist_images, mnist_labels), _ = tf.keras.datasets.mnist.load_data()
mnist_images, mnist_labels = mnist_images[:128], mnist_labels[:128]
dataset = tf.data.Dataset.from_tensor_slices((tf.cast(mnist_images[..., tf.newaxis] / 255, tf.float32), tf.cast(mnist_labels, tf.int64)))
dataset = dataset.shuffle(1000).batch(32)

tot_loss = []
epochs = 1

mnist_model = tf.keras.Sequential([
    tf.keras.layers.Conv2D(16, [3, 3], activation='relu'),
    tf.keras.layers.Conv2D(16, [3, 3], activation='relu'),
    tf.keras.layers.GlobalAveragePooling2D(),
    tf.keras.layers.Dense(10)
])

optimizer = tf.keras.optimizers.Adam()
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

for epoch in range(epochs):
    loss_history = []
    for (batch, (images, labels)) in enumerate(dataset):
        with tf.GradientTape() as tape:
            logits = mnist_model(images, training=True)
            loss_value = loss_object(labels, logits)
        loss_history.append(loss_value.numpy().mean())
        grads = tape.gradient(loss_value, mnist_model.trainable_variables)
        optimizer.apply_gradients(zip(grads, mnist_model.trainable_variables))
    tot_loss.append(np.mean(loss_history))

tot_loss
